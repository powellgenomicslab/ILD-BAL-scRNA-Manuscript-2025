setwd("/directflow/SCCGGroupShare/projects/petall/ild-pathwayAnalysis-bal-scRNA/")

library(tidyverse)
library(scPagwas)

args <- commandArgs(trailingOnly = TRUE)

pathway <- args[1]

gwas_path <- "data/01-raw/allen2022_metaGwas5way_sumstats.txt"

if(pathway == "reactome"){
    output_dir <- "data/03-final/scpagwas_reactome_allen2022_ipf5way"

  # When running subsequent analyses (e.g. same scRNA data, different gwas data) with Pagwas (you can specify this and use the same single cell pathway analysis step)
    pagwas <- scPagwas_main(Pagwas = NULL,
                            gwas_data = gwas_path,
                            Single_data = "data/02-preprocessed/BAL_integrated_scPagwas.rds",
                            output.prefix = "IPF",
                            output.dirs = output_dir,
                            n.cores = 6,
                            Pathway_list = genes.by.reactome.pathway, #Genes_by_pathway_kegg,
                            assay = "RNA",
                            singlecell = TRUE,
                            iters_singlecell = 100,
                            celltype = TRUE,
                            block_annotation = block_annotation,
                            chrom_ld = chrom_ld)
} else {
    output_dir <- "data/03-final/scpagwas_kegg_allen2022_ipf5way"

    pagwas <- scPagwas_main(Pagwas = NULL,
                            gwas_data = gwas_path,
                            Single_data = "data/02-preprocessed/BAL_integrated_scPagwas.rds",
                            output.prefix = "allen2022_IPF",
                            output.dirs = output_dir,
                            n.cores = 6,
                            Pathway_list = Genes_by_pathway_kegg,
                            assay = "RNA",
                            singlecell = TRUE,
                            iters_singlecell = 100,
                            celltype = TRUE,
                            block_annotation = block_annotation,
                            chrom_ld = chrom_ld)
}


saveRDS(pagwas, file = paste0(output_dir, "/pagwas.rds"))


# Running for 5 way
# gwas_path <- "/g/data/ei56/pa3687/ILD_scLinker/data/plink-output/IPF_EUR_5way2020_prune_gwas_data.txt"
#
# output_dir <- paste0("scpagwas_", str_split(gwas_path, "/")[[1]][9]) %>%
#   gsub("_prune_gwas_data.txt", "", .)
#
# pagwas_5way <- scPagwas_main(Pagwas = pagwas_3way,
#                         gwas_data = gwas_path,
# #                        Single_data = "/g/data/ei56/pa3687/ILD_scLinker/output/BAL_integrated_scPagwas.rds",
#                         output.prefix = "IPF",
#                         output.dirs = output_dir,
#                         n.cores = 6,
#                         Pathway_list = Genes_by_pathway_kegg,
#                         assay = "RNA",
#                         singlecell = TRUE,
#                         iters_singlecell = 100,
#                         celltype = TRUE,
#                         block_annotation = block_annotation,
#                         chrom_ld = chrom_ld)

# > pagwas_5way <- scPagwas_main(Pagwas = pagwas_3way,
# +                         gwas_data = gwas_path,
# + #                        Single_data = "/g/data/ei56/pa3687/ILD_scLinker/output/BAL_integrated_scPagwas.rds",
# +                         output.prefix = "IPF",
# +                         output.dirs = output_dir,
#            n.cores +                         n.cores = 6,
# +                         Pathway_list = Genes_by_pathway_kegg,
# +                         assay = "RNA",
# +                         singlecell = TRUE,
# +                         iters_singlecell = 100,
# +                         celltype = TRUE,
# +                         block_annotation = block_annotation,
# +                         chrom_ld = chrom_ld)
# ##------ Mon Sep  2 09:11:02 2024 ------## ******* 1st: Single_data_input function start! ********
# done!
# ##------ Mon Sep  2 09:27:15 2024 ------## ******* 3rd: GWAS_summary_input function start! ********
# ** Start to read the gwas_data!
# |--------------------------------------------------|
# |==================================================|
# Input gwas summary data frame!
# No "chr" from chrom!, now pasting it!
# Filtering out SNPs with MAF criterion
# done!
# ##------ Mon Sep  2 09:27:40 2024 ------## ******* 4th: SnpToGene start!! ********
# ##------ Mon Sep  2 09:27:44 2024 ------## ******* 5th: Pathway_annotation_input function
#                 start! ********
# Error in Pathway_annotation_input(Pagwas = Pagwas, block_annotation = block_annotation) :
#   no match for Pathway gene and VariableFeatures
# >
