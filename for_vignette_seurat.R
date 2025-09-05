library( ggplot2 )
library( Seurat )
library( exprest )

seu <- LoadSeuratRds("Blood_TSP1_30_version2d_10X_smartseq_scvi_Nov122024_SeuratObject.rds")

seu <- NormalizeData(seu)
VlnPlot(seu, "LCP1", group.by="cell_ontology_class", alpha=.1 ) + guides(fill = "none")

a <- est_dens_for_groups( LayerData(seu,"counts")["LCP1",], seu$nCount_RNA, seu$cell_ontology_class )

ggplot(a) + geom_violin(aes( x=group, y=log10lambda, weight=density, fill=group )) +
  guides(fill="none") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot( seu@meta.data[,c("cell_ontology_class","nCount_RNA")] ) +
  geom_violin(aes(x=cell_ontology_class,y=log10(nCount_RNA),fill=cell_ontology_class)) +
  guides(fill="none") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
