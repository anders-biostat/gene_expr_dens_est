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

FindClusters(seu,resolution=0.1)-> seu

#Rcpp::sourceCpp("simplex_sampler_threaded.cc", verbose=TRUE, rebuild=TRUE ); 

Sys.setenv(PKG_CXXFLAGS = system("pkg-config --cflags libcmaes", intern=TRUE) )
Sys.setenv(PKG_LIBS = system("pkg-config --libs libcmaes", intern=TRUE) )
Rcpp::sourceCpp("test_cmaes.cc", verbose=TRUE, rebuild=TRUE)

estimate_density_mcmc <- function( k, s ) {
   
   # prepare gamma grid
   stepsize <- log(1.125)  # log(nu)
   mu <- exp( seq( log(3e-5), 0, by=stepsize ) )
   shape <- 10 / stepsize  # kappa
   scale <- mu / shape    # theta
   
   # pre-compute likelihoods
   pmf_nb <- sapply( mu, function(mu_) dnbinom( k, mu=s*mu_, size=shape ) )
   
   # run simulated annealing
   draws <- sample_mixture_weights_threads( pmf_nb, n_burnin_draws=0L, n_keep_draws=10000L, 
         temperature_decrease=.95, global_seed = as.integer( runif(1,0,100000) ),  )
   
   # Get mixture weights
   spi <- sapply( draws, function(x) x[ nrow(x), ] )
   spim <- rowMeans(spi)
   
   ans <- list( log_lambda =  seq( -6, 0, length.out=1000 ) )
   ans$density <- sapply( 1:length(mu), function(i)
      dgamma( 10^ans$log_lambda, shape=mu[i]/scale[i], scale=scale[i] )*(10^ans$log_lambda)*log(10) ) %*% spim
   ans$error_est <- max( dist(t(spi)) )

   ans
}

estimate_density <- function( k, s ) {
   
   # prepare gamma grid
   stepsize <- log(1.125)  # log(nu)
   mu <- exp( seq( log(3e-5), 0, by=stepsize ) )
   shape <- 10 / stepsize  # kappa
   scale <- mu / shape    # theta
   
   # pre-compute likelihoods
   pmf_nb <- sapply( mu, function(mu_) dnbinom( k, mu=s*mu_, size=shape ) )

   x <- do_sann( pmf_nb )
   
   # Get mixture weights
   v <- 1/(1+exp(-x))
   spi <- c( 1, cumprod(v) ) * c( 1-v, 1 ) 

   dens_est <- data.frame(
      log10lambda = seq( -6, 0, length.out=1000 ) )
   dens_est$density <- (
      sapply( 1:length(mu), function(i)
         dgamma( 10^dens_est$log10lambda, shape=mu[i]/scale[i], scale=scale[i] ) * 
            (10^dens_est$log10lambda)*log(10) ) 
      %*% spi )

   dens_est
}


gene <- "CD8A"

VlnPlot(seu,gene)


lapply( 
   levels( seu$seurat_clusters ),
   function(cl) {
      cat(cl)
      estimate_density( 
         LayerData(seu,"counts")[ gene, seu$seurat_clusters==cl ],
         colSums( LayerData(seu,"counts") )[seu$seurat_clusters==cl] ) } ) -> ans
names(ans) <- levels( seu$seurat_clusters )
max( sapply( ans, function(a) a$error_est) )

plot( NULL, ylim=c(-5.3,-1), xlim=c(0,7), main=gene, ylab="log10 expression fraction", xlab="Cluster" )
for( cl in names(ans) ) {
   x <- ans[[cl]]$log10lambda
   y <- ans[[cl]]$density
   y <- ifelse( y > sum(y)*1e-3, y, 0 )
   polygon( 
      as.numeric(cl) + c( y, -rev(y) )/15,
      c( x, rev(x) ),
      border=NA, col="#ffb0b0c0" )
}
   
for( cl in names(ans) ) {

   k <- LayerData(seu,"counts")[ gene, seu$seurat_clusters==cl ]
   s <- colSums( LayerData(seu,"counts") )[seu$seurat_clusters==cl]
   mu <- mean(k/s)
   v <- var(k/s) - mu * mean(1/s)

   x <- ans[[cl]]$log_lambda
   y <- dgamma( 10^ans[[cl]]$log_lambda, shape = mu^2/v, scale = v/mu ) * (10^ans[[cl]]$log_lambda)*log(10)
   y <- ifelse( y > sum(y)*1e-3, y, 0 )
   polygon( 
      as.numeric(cl) + c( y, -rev(y) )/5,
      c( x, rev(x) ),
      border=NA, col="#00FF5010" )
}



cl <- "1"
k <- LayerData(seu,"counts")[ gene, seu$seurat_clusters==cl ]
s <- colSums( LayerData(seu,"counts") )[seu$seurat_clusters==cl]
mu <- mean(k/s)
v <- var(k/s) - mu * mean(1/s)

plot( ans[[cl]]$log_lambda, ans[[cl]]$density, type="l" )
x <- seq(-20,0,length.out=1000)
lines( x, dgamma( 10^x, shape = mu^2/v, scale = v/mu ) * (10^x)*log(10), col="red" )

