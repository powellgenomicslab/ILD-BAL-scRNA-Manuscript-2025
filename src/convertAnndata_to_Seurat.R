library(Seurat)
library(reticulate)
library(sceasy)


use_condaenv('ild-bal-scrna')
h5ad_file <- 'data/03-final/integrated.h5ad'
sceasy::convertFormat(h5ad_file, from="anndata", to="seurat",
                       outFile='data/03-final/integrated.RDS')


bal_ild <- readRDS("data/03-final/integrated.RDS")
bal_ild

# Ensure UMAP looks the same between python/R --  confirmed
png("reports/figures/dimplot_seurat.png", height=20, width=35, units="in", res=200)
DimPlot(bal_ild, reduction = "umap", group.by = c("Sample", "disease_final"))
dev.off()

# Clustering using Seurat on the scVI embedding
bal_ild <- FindNeighbors(bal_ild, reduction="scVI", dims=1:10)
bal_ild <- FindClusters(bal_ild, resolution=2, cluster.name = "scvi_clusters")


png("reports/figures/dimplot_scvi_clusters_res2.png", height=10, width=15, units="in", res=200)
DimPlot(bal_ild, reduction = "umap", group.by = c("disease_final", "scvi_clusters"))
dev.off()


png("reports/figures/featurePlot_cd4_cd8.png", height=10, width=15, units="in", res=200)
FeaturePlot(bal_ild, features = c("CD4", "CD8A", "CD8B"), max.cutoff=5)
dev.off()
