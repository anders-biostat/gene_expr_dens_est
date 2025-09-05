# This script converts a h5ad objects as downloaded, eg., from the Tabula Sapiens FigShare repo, to a Seurat object
#    Tabula Sapiens repo on FigShare: https://figshare.com/articles/dataset/Tabula_Sapiens_v2/27921984

library(Matrix)
library(tidyverse)
library(reticulate)
use_python("/usr/bin/python")
library(Seurat)
anndata <- import("anndata")

ad = anndata$read_h5ad( "../Blood_TSP1_30_version2d_10X_smartseq_scvi_Nov122024.h5ad" )

# Convert sparse integer matrix (auto-conversion does not work)
a <- py_eval("r['ad'].layers['raw_counts']", convert = FALSE)
raw_counts <- sparseMatrix( i=py_to_r(a$indices), p=py_to_r(a$indptr), x=py_to_r(a$data), index1=FALSE )

# Set dimnames
colnames(raw_counts) <- rownames(ad$obs)
rownames(raw_counts) <- rownames(ad$var)

# Make Seurat object, include cell metadata
seu <- CreateSeuratObject( raw_counts, meta.data=ad$obs, project="TabulaSapiensBlood" )

# Copy gene metadata
seu[["RNA"]]@meta.data <- ad$var

# Copy UMAP result
ump <- ad$obsm[["X_umap"]]
colnames(ump) <- paste0( "UMAP_", 1:2 )
rownames(ump) <- colnames(seu)
seu[["umap"]] <- CreateDimReducObject( embeddings=ump, assay="RNA", key="UMAP_", global=TRUE )

# Save object
SaveSeuratRds( seu, "Blood_TSP1_30_version2d_10X_smartseq_scvi_Nov122024_SeuratObject.rds" )
