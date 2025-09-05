library(Matrix)
library(tidyverse)
library(reticulate)
use_python("/usr/bin/python")
library(Seurat)

anndata <- import("anndata")
ad = anndata$read_h5ad( "../Blood_TSP1_30_version2d_10X_smartseq_scvi_Nov122024.h5ad" )
ad$layers[["raw_counts"]]$data
a <- py_eval("r['ad'].layers['raw_counts']", convert = FALSE)
raw_counts <- sparseMatrix( i=py_to_r(a$indices), p=py_to_r(a$indptr), x=py_to_r(a$data), index1=FALSE )
colnames(raw_counts) <- rownames(ad$obs)
rownames(raw_counts) <- rownames(ad$var)
seu <- CreateSeuratObject( raw_counts, meta.data=ad$obs, project="TabulaSapiensBlood" )
seu[["RNA"]]@meta.data <- ad$var

ump <- ad$obsm[["X_umap"]]
colnames(ump) <- paste0( "UMAP_", 1:2 )
rownames(ump) <- colnames(seu)
seu[["umap"]] <- CreateDimReducObject( embeddings=ump, assay="RNA", key="UMAP_", global=TRUE )

seu <- NormalizeData(seu)
VlnPlot(seu, "LCP1", group.by="cell_ontology_class", alpha=.1 ) + guides(fill = "none")

SaveSeuratRds( seu, )

library( exprest )


stopifnot( all( seu$nCount_RNA == colSums(LayerData(seu,"counts")) ) )

k <- LayerData(seu,"counts")["LCP1",]
s <- seu$nCount_RNA
g <- seu$cell_ontology_class
pb <- txtProgressBar(max=length(s), style=3)
a <- tapply( 1:length(s), g, function(cells) {
  if( length(cells)>1 )
    a <- exprest::est_dens( k[cells], s[cells] ) 
  else
    a <- NA
  setTxtProgressBar( pb, getTxtProgressBar(pb) + length(cells) ) 
  a }, simplify=FALSE )

do.call( bind_cols, a )

a2 %>% map( as_tibble ) %>% bind_rows(.id="group") %>%
mutate(density=drop(density)) %>%
ggplot() + geom_violin(aes( x=group, y=log10lambda, weight=density, fill=group )) +
guides(fill="none") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

seu@meta.data[,c("cell_ontology_class","nCount_RNA")] %>%
ggplot() + geom_violin(aes(x=cell_ontology_class,y=log10(nCount_RNA),fill=cell_ontology_class)) +
guides(fill="none") + theme(axis.text.x = element_text(angle = 45, hjust = 1))



k <- LayerData(seu,"counts")["LCP1",]
s <- seu$nCount_RNA
g <- seu$cell_ontology_class
max_per_group <- 700
pb <- txtProgressBar(max=sum( pmin( table(g), max_per_group ) ), style=3)
system.time( a2 <- tapply( 1:length(s), g, function(cells) {
  if( length(cells)>1 ) {
    a <- est_dens( k[cells], s[cells] ) 
  }
  else
    a <- NA
  setTxtProgressBar( pb, getTxtProgressBar(pb) + min( 700, length(cells) ) )
  a }, simplify=FALSE ) )
