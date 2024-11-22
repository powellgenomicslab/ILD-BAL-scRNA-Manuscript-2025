library(Seurat)
library(tidyverse)


# print("Convert from Scanpy to Seurat...")
# bal_ild <- read_h5ad("output/BAL_integrated_annotated.h5ad")

# umap_coord <- bal_ild$obsm$X_umap
# rownames(umap_coord) <- rownames(bal_ild)
# umap_coord <- as.matrix(umap_coord)
# colnames(umap_coord) <- c("UMAP_1", "UMAP_2")

# bal_ild <- CreateSeuratObject(counts = t(bal_ild$X), meta.data = bal_ild$obs)
# pagwas_3way[['umap']] <- CreateDimReducObject(embeddings = umap_coord, key = "UMAP_", global = T, assay = "RNA")
BAL_allen_normalised_reannotated <- readRDS("data/02-preprocessed/BAL_allen_normalised_reannotated.rds")
cell_meta <- read_csv("data/02-preprocessed/corrected_metadata.csv")
# macrophage_meta <- read_csv("data/02-preprocessed/mdm_cluster_metadata.csv")

cell_meta$IPF[match(rownames(BAL_allen_normalised_reannotated@meta.data), cell_meta$...1)] -> BAL_allen_normalised_reannotated@meta.data$IPF
cell_meta$telomere_length[match(rownames(BAL_allen_normalised_reannotated@meta.data), cell_meta$...1)] -> BAL_allen_normalised_reannotated@meta.data$telomere_length
cell_meta$disease_final[match(rownames(BAL_allen_normalised_reannotated@meta.data), cell_meta$...1)] -> BAL_allen_normalised_reannotated@meta.data$disease_final

# macrophage_meta$mdm_clusters <- paste0("mdm_", macrophage_meta$mdm_clusters)
# matching_barcodes <- na.omit(match(macrophage_meta$...1,rownames(BAL_allen_normalised_reannotated@meta.data)))
# BAL_allen_normalised_reannotated@meta.data$Cell_Subtype <- as.character(BAL_allen_normalised_reannotated@meta.data$Cell_Subtype)
# BAL_allen_normalised_reannotated@meta.data$Cell_Subtype[matching_barcodes] <- macrophage_meta$mdm_clusters[na.omit(match(rownames(BAL_allen_normalised_reannotated@meta.data), macrophage_meta$...1))]

Idents(object = BAL_allen_normalised_reannotated) <- "Cell_Subtype"
BAL_allen_normalised_reannotated <- NormalizeData(BAL_allen_normalised_reannotated, normalization.method = "LogNormalize", scale.factor = 10000)
BAL_allen_normalised_reannotated <- ScaleData(BAL_allen_normalised_reannotated)
saveRDS(BAL_allen_normalised_reannotated, "data/02-preprocessed/BAL_integrated_scPagwas.rds")


# ## Split into another script post scPagwas
# bal_ild <- readRDS("output/BAL_integrated_scPagwas.rds")
# bal_ild_meta <- bal_ild@meta.data

# ## read in scPagwas results
# ipf_3way_2019_scpagwas <- read_csv("/g/data/ei56/pa3687/ILD_scLinker/IPF_EUR_3way2019/IPF_singlecell_scPagwas_score_pvalue.Result.csv")

# all(rownames(bal_ild@meta.data) == ipf_3way_2019_scpagwas[,1])

# bal_ild$scPagwas.TRS.Score <- ipf_3way_2019_scpagwas$scPagwas.TRS.Score
# bal_ild$scPagwas.downTRS.Score <- ipf_3way_2019_scpagwas$scPagwas.downTRS.Score
# bal_ild$scPagwas.gPAS.score <- ipf_3way_2019_scpagwas$scPagwas.gPAS.score
# bal_ild$Random_Correct_BG_p <- ipf_3way_2019_scpagwas$Random_Correct_BG_p
# bal_ild$Random_Correct_BG_adjp <- ipf_3way_2019_scpagwas$Random_Correct_BG_adjp
# bal_ild$Random_Correct_BG_z <- ipf_3way_2019_scpagwas$Random_Correct_BG_z
# bal_ild$Random_Correct_BG_adjp_0.05 <- ifelse(bal_ild$Random_Correct_BG_adjp<0.05, TRUE, FALSE)

# png("umap-bal_ild.png")
# DimPlot(bal_ild)
# dev.off()


# png("umap-bal_ild-TRS.png")
# FeaturePlot(object = bal_ild, features = 'scPagwas.TRS.Score')
# dev.off()


# png("umap-bal_ild-trs-sig.png")
# FeaturePlot(object = bal_ild, features = 'Random_Correct_BG_adjp_0.05')
# dev.off()
