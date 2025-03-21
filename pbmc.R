library( tidyverse )
library( Seurat )

Read10X_h5( "~/Downloads/10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.h5" ) %>%
CreateSeuratObject() %>%
NormalizeData() %>%
FindVariableFeatures() %>%
ScaleData() %>%
RunPCA() %>%
FindNeighbors() %>%
FindClusters() %>%
RunUMAP(dims=1:50) -> seu
UMAPPlot(seu, label=TRUE) + coord_equal()

VlnPlot(seu,"IL7R")

Rcpp::sourceCpp("simplex_sampler_threaded.cc", verbose=TRUE, rebuild=TRUE ); 

estimate_density <- function( k, s ) {
   
   # preparegamma grid
   stepsize <- log(1.05)  # log(nu)
   mu <- exp( seq( log(3e-5), 0, by=stepsize ) )
   shape <- 30 / stepsize  # kappa
   scale <- mu / shape    # theta
   
   # pre-compute likelihoods
   pmf_nb <- sapply( mu, function(mu_) dnbinom( k, mu=s*mu_, size=shape ) )
   
   # run simulated annealing
   draws <- sample_mixture_weights_threads( pmf_nb, n_burnin_draws=2000L, n_keep_draws=100000L, 
         temperature_decrease=1.0, global_seed = as.integer( runif(1,0,100000) ),  )
   
   # Get mixture weights
   spi <- sapply( draws, function(x) x[ nrow(x), ] )
   spim <- rowMeans(spi)
   
   ans <- list( log_lambda =  seq( -6, 0, length.out=1000 ) )
   ans$density <- sapply( 1:length(mu), function(i)
      dgamma( 10^ans$log_lambda, shape=mu[i]/scale[i], scale=scale[i] )*(10^ans$log_lambda)*log(10) ) %*% spim
   ans$error_est <- max( dist(t(spi)) )

   ans
}

gene <- "IL7R"
lapply( 
   levels( seu$seurat_clusters ),
   function(cl) {
      cat(cl)
      estimate_density( 
         LayerData(seu,"counts")[ gene, seu$seurat_clusters==cl ],
         colSums( LayerData(seu,"counts") )[seu$seurat_clusters==cl] ) } ) -> ans
names(ans) <- levels( seu$seurat_clusters )
max( sapply( ans, function(a) a$error_est) )

plot( NULL, ylim=c(-5,-1), xlim=c(0,21) )
for( cl in names(ans) ) {
   x <- ans[[cl]]$log_lambda
   y <- ans[[cl]]$density
   y <- ifelse( y > sum(y)*1e-3, y, 0 )
   polygon( 
      as.numeric(cl) + c( y, -rev(y) )/7,
      c( x, rev(x) ),
      border=NA, col="#b0b0b080" )
}
   
