library(tidyverse)

# Read and preprocess the GWAS data
gwas <- read_delim("/g/data/ei56/pa3687/IPF_Bothsex_inv_var_meta_GBMI_with_Allen_021121.txt.gz")

gwas <- gwas %>%
  rename("chrom" = "#CHR",
         "pos" = "POS",
         "beta" = "inv_var_meta_beta",
         "se" = "inv_var_meta_sebeta",
         "maf" = "all_meta_AF") %>%
  select(chrom, pos, REF, ALT, rsid, beta, se, maf)

# ---

# Write processed GWAS data to a file
write.table(gwas, file = "data/IPF_bothsex_gbmi-allen2020.txt", row.names = FALSE, quote = FALSE)

# ---

# Extract SNP list from the GWAS file
system("awk '{print $5 }' data/IPF_bothsex_gbmi-allen2020.txt > data/IPF_bothsex_gbmi-allen2020_SNP_list.txt")

# Loop over chromosomes 1 to 22 to run plink commands
for (i in 1:22) {
    cat("chr:", i, "\n")
    system(sprintf("plink --bfile code/scLinker/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.%d --extract data/IPF_bothsex_gbmi-allen2020_SNP_list.txt --make-bed --out data/plink-output/1000G.EUR.QC.IPF_bothsex_gbmi-allen2020_%d_filtered", i, i))
}

# Loop over chromosomes 1 to 22 to perform LD pruning with plink
for (i in 1:22) {
    cat("chr:", i, "\n")
    system(sprintf("plink --bfile data/plink-output/1000G.EUR.QC.IPF_bothsex_gbmi-allen2020_%d_filtered --indep-pairwise 50 5 0.8 --out data/plink-output/IPF_bothsex_gbmi-allen2020_%d_plink_prune_filtered_LD0.8", i, i))
}

# Concatenate prune files into one
system("cat data/plink-output/IPF_bothsex_gbmi-allen2020_*.prune.in > data/plink-output/IPF_bothsex_gbmi-allen2020_LD0.8.prune")

# Read the pruned SNP list and join with the original GWAS data
gwas <- read_delim("data/IPF_bothsex_gbmi-allen2020.txt")
SNP_prune <- read_table("data/plink-output/IPF_bothsex_gbmi-allen2020_LD0.8.prune", col_names=F)
SNP_prune <- SNP_prune[!duplicated(unlist(SNP_prune)),]
colnames(SNP_prune) <- "rsid"

# Perform inner join
gwas <- gwas %>% inner_join(SNP_prune, by="rsid")
print(nrow(gwas))

# Write the filtered GWAS data to a file
write.table(gwas, file="data/plink-output/IPF_bothsex_gbmi-allen2020_prune_gwas_data.txt", row.names=F, quote=F)
