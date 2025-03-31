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


library( exprest )

gene <- "CD8A"

VlnPlot(seu,gene)


lapply( 
   levels( seu$seurat_clusters ),
   function(cl) {
      cat(cl)
      est_dens( 
         LayerData(seu,"counts")[ gene, seu$seurat_clusters==cl ],
         colSums( LayerData(seu,"counts") )[seu$seurat_clusters==cl] ) } ) -> ans
names(ans) <- levels( seu$seurat_clusters )
#max( sapply( ans, function(a) a$error_est) )

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

   x <- ans[[cl]]$log10lambda
   y <- dgamma( 10^ans[[cl]]$log10lambda, shape = mu^2/v, scale = v/mu ) * (10^ans[[cl]]$log10lambda)*log(10)
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

