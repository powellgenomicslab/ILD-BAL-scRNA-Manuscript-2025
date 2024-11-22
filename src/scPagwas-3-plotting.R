library(tidyverse)
library(scPagwas)
library(data.table)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
outDir <- args[1]


pagwas_path <- paste0(outDir, "/pagwas.rds")

pagwas <- readRDS(pagwas_path)
## pagwas <- readRDS("/g/data/ei56/pa3687/ILD_scLinker/output/scpagwas_IPF_EUR_5way2020.rds")


pval_threshold <- quantile(pagwas$Random_Correct_BG_adjp, 0.25)

scPagwas_Visualization(Single_data=pagwas,
                      p_thre = pval_threshold,
                      FigureType = "umap",
                      width = 7,
                      height = 7,
                      lowColor = "white", 
                      highColor = "red",
                      output.dirs=outDir,
                      size = 0.5,
                      do_plot = F)

# Proportion Plot of Positive vs Negative Cells

## Add meta
pagwas$positiveCells <- rep(0, ncol(pagwas))
pagwas$positiveCells[pagwas$Random_Correct_BG_adjp < pval_threshold] <- 1

## Identify numbers for each cell type
process_data <- function(meta_data, ident_column, group_column) {
  # Create data table from meta data
  dt <- data.table(
    "ident" = as.character(meta_data[[ident_column]]),
    "group" = as.character(meta_data[[group_column]])
  )
  
  # Create n_ident column
  dt[, n_ident := paste0(ident, " (n=", .N, ")"), by = ident]
  
  # Generate factor levels based on numeric ordering
  vec_factorLevels <- dt$n_ident[gsub("\\ .*", "", dt$n_ident) %>%
                                    as.numeric() %>%
                                    order()] %>%
    unique()
  
  # Make n_ident a factor with ordered levels
  dt[, n_ident := factor(n_ident, levels = vec_factorLevels, ordered = TRUE)]
  
  # Calculate summary counts
  dt_sum <- dt[, .N, by = .(n_ident, group)]
  
  return(dt_sum)
}

# Function for filtering by "IPF" and applying process_data
process_ipf_data <- function(meta_data) {
  # Filter for IPF
  meta_data_ipf <- meta_data[meta_data$IPF == "IPF", ]
  
  # Apply process_data on filtered data (assuming 'Cell_Subtype' and 'positiveCells' columns exist)
  dt_sum_ipf <- process_data(meta_data_ipf, "Cell_Subtype", "positiveCells")
  
  return(dt_sum_ipf)
}

# Example of applying functions
dt_sum_all <- process_data(pagwas@meta.data, "Cell_Subtype", "positiveCells")
dt_sum_ipf <- process_ipf_data(pagwas@meta.data)

## ggplot for all data
p_all <- ggplot(
  dt_sum_all,
  aes(x = n_ident, y = N, fill = factor(group))
  ) +
  geom_bar(
    position = "fill",
    stat = "identity",
    width = 0.6,
    show.legend = FALSE
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("0" = "lightgray",
                 "1" = "blue")) +
  theme_bw() +
  theme(
    axis.title.x = element_text(vjust = 0),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.margin = margin(t=0.1, b=0.1, r=0.1, l=2, "cm")
  ) +
  labs(
    x = "Cell type", 
    y = "Proportion", 
    fill = "Positive Cells",
    title = "Proportion of Positive vs Negative Cells (All Data)"
  )

ggsave(plot = p_all, filename = paste0(outDir, "/cellProportionPathway_all.png"), width = 7, height = 7, dpi = 300)

## ggplot for IPF data
p_ipf <- ggplot(
  dt_sum_ipf,
  aes(x = n_ident, y = N, fill = factor(group))
  ) +
  geom_bar(
    position = "fill",
    stat = "identity",
    width = 0.6,
    show.legend = FALSE
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("0" = "lightgray",
                 "1" = "blue")) +
  theme_bw() +
  theme(
    axis.title.x = element_text(vjust = 0),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.margin = margin(t=0.1, b=0.1, r=0.1, l=2, "cm")
  ) +
  labs(
    x = "Cell type", 
    y = "Proportion", 
    fill = "Positive Cells",
    title = "Proportion of Positive vs Negative Cells (IPF)"
  )

ggsave(plot = p_ipf, filename = paste0(outDir, "/cellProportionPathway_ipf.png"), width = 7, height = 7, dpi = 300)
# ---

# Correlation Scatterplot for Pearson Correlation Coefficients (PCC)

cor_df <- data.frame(genes = rownames(pagwas@misc$PCC), PCC = pagwas@misc$PCC[, 1], text = rep(NA, nrow(pagwas@misc$PCC)))

cor_df <- cor_df[order(cor_df$PCC, decreasing = T), ]
cor_df$text[1:10] <- cor_df$genes[1:10]

cor_df$order <- seq_len(nrow(cor_df))

p <- ggplot(data = cor_df) +
  geom_point(mapping = aes(x = order, y = PCC, color = PCC)) +
  scale_colour_gradient2(
    low = "#035397",
    mid = "white",
    high = "#F32424",
    midpoint = 0
  ) +
  theme_classic() +
  ggrepel::geom_text_repel(aes(
    x = order, y = PCC,
    label = text
  ),
  max.overlaps = 20,
  na.rm = T,
  force = 8,
  force_pull = 0.1,
  )

ggsave(plot = p, filename = paste0(outDir, "/correlation_scatter.png"), width = 7, height = 7, dpi = 300)

# ---

# Cell-type bootstrap Plots
## Bootstrap_P_Barplot(p_results=pagwas_3way@misc$bootstrap_results$bp_value[-1],
##   p_names=rownames(pagwas_3way@misc$bootstrap_results)[-1],
##   figurenames = "Bootstrap_P_Barplot.pdf",
##   width = 5,
##   height = 7,
##   do_plot=T,
##   title = "ipf3way_bal")
##
## Bootstrap_estimate_Plot(bootstrap_results=pagwas_3way@misc$bootstrap_results,
## ,        figurenames = "estimateplot.pdf",
## ,        width = 9,
## ,        height = 7,
## ,        do_plot=T)
# ---

# Pathway plot

## Proportion
proportion_list <- tapply(
as.vector(Idents(pagwas)),
Idents(pagwas), function(x) {
    scPagwasPaHeritability <- t(GetAssayData(pagwas,
                                            assay = "scPagwasPaHeritability"))
    a <- apply(scPagwasPaHeritability, 2,
                function(y) sum(y > 0) / length(y))
    return(unlist(a))
}
)

proportion_df <- Reduce(function(dtf1, dtf2) cbind(dtf1, dtf2),
                        proportion_list)
colnames(proportion_df) <- names(proportion_list)

## Rank p-value
scPathrankP <- -log10(pagwas@misc$scPathways_rankPvalue + 1e-20)

top_function <- function(para_mat, n_path_to_keep) {
para_mat$path <- rownames(para_mat)

## Only keep genes with a unique name and tidy data.
para_mat <- para_mat %>%
    add_count(path) %>%
    filter(n == 1) %>%
    select(-n) %>%
    gather(key = celltypes, value = paras, -path) %>%
    as_tibble()

para_mat <- para_mat %>%
    group_by(celltypes) %>%
    mutate(paras_sum_mean = (paras * nrow(scPathrankP)) / sum(paras))

para_mat <- para_mat %>%
    group_by(path) %>%
    mutate(specificity = (paras_sum_mean * nrow(scPathrankP)) / sum(paras_sum_mean)) %>%
    ungroup()

d_spe <- para_mat %>% filter(paras_sum_mean > 1)
d_spe <- d_spe %>%
    group_by(celltypes) %>%
    top_n(., n_path_to_keep, specificity)
return(d_spe)
}

spe2 <- top_function(para_mat = scPathrankP,
                     n_path_to_keep = 5)

celltypes <- levels(as.factor(pagwas$Cell_Subtype))
spe2 <- spe2[spe2$celltypes %in% celltypes, ]
spe <- unique(spe2$path)

## Merge to plot
rownames(proportion_df) <- gsub("-", "_", rownames(proportion_df))

spe<-intersect(rownames(proportion_df),spe)
scPathrankP <- scPathrankP[spe, celltypes]
proportion_df <- proportion_df[spe, celltypes]
proportion_df <- as.data.frame(proportion_df)

scPathrankP[scPathrankP > 15] <- 15
scPathrankP[scPathrankP < -log10(0.05)] <- 0

proportion_df$pathways <- rownames(proportion_df)
scPathrankP$pathways <- rownames(scPathrankP)

gg_proportion <- reshape2::melt(proportion_df,
                                id.vars = "pathways",
                                variable.name = "celltypes",
                                value.name = "proportion")

gg_rankp <- reshape2::melt(scPathrankP,
                            id.vars = "pathways",
                            variable.name = "celltypes",
                            value.name = "logrankvalue")

gg_dot <- merge(gg_rankp, gg_proportion)

## Plot the pathway dot plot
library(ggdendro)

pathway_wide <- gg_dot
pathway_wide <- pathway_wide %>%
  pivot_wider(
    names_from = celltypes,
    values_from = c(logrankvalue, proportion)
  )

rownames(pathway_wide) <- pathway_wide$pathways
pathway_dendro <- as.dendrogram(hclust(d = dist(x = pathway_wide)))


p <- ggdendrogram(data = pathway_dendro, rotate = TRUE) +
  theme(axis.text.y = element_text(size = 6))
ggsave(plot = p, filename = paste0(outDir,"/pathway_dendrogram.png"), height = 12, width = 8, dpi = 300)

### Extract the order of the tips in the dendrogram
pathway_order <- order.dendrogram(pathway_dendro)
### Order the levels according to their position in the cluster
gg_dot$pathways <- factor(x = gg_dot$pathways,
                               levels = pathway_wide$pathways[pathway_order],
                               ordered = TRUE)

if(any(grep("REACTOME", gg_dot$pathways))){
   dotplot <- ggplot(gg_dot, aes(x = celltypes, y = pathways)) +
    geom_point(data = subset(gg_dot, logrankvalue != 0),
               aes(size = logrankvalue, color = proportion)) +
    scale_size_continuous(name = "Log Rank Value") +
    scale_color_gradient(name = "Proportion", low = "grey", high = "brown") +
    labs(x = "Cell Types", y = "Pathways") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          plot.margin = margin(t=0.8, b=0.8, r=0.5, l=2, "cm"))

    ggsave(plot = dotplot, filename = paste0(outDir, "/pathway_plot.png"), height = 12, width = 16, dpi = 300)
} else {
    dotplot <- ggplot(gg_dot, aes(x = celltypes, y = pathways)) +
    geom_point(data = subset(gg_dot, logrankvalue != 0),
               aes(size = logrankvalue, color = proportion)) +
    scale_size_continuous(name = "Log Rank Value") +
    scale_color_gradient(name = "Proportion", low = "grey", high = "brown") +
    labs(x = "Cell Types", y = "Pathways") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")

    p <- p + theme(plot.margin = margin(t=1, b=2.8, r=0.1, l=0.1, "cm"))

    ## Combine the plots
    combined_plot <- cowplot::plot_grid(dotplot, p, ncol = 2, align = "v", rel_heights = c(1, 1))

    ggsave(plot = combined_plot, filename = paste0(outDir, "/pathway_plot.png"), height = 12, width = 16, dpi = 300)

}

# ---

# Plot TRS Across Cell Types
# Violin plot of TRS Score by Cell Subtype
# Violin plot of TRS Score by Cell Subtype
p_trs <- ggplot(pagwas@meta.data, aes(x = Cell_Subtype, y = scPagwas.downTRS.Score2, color = Cell_Subtype)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t=0.1, b=0.1, r=0.1, l=2, "cm")
  ) +
  labs(
    x = "Cell Type",
    y = "TRS Score",
    title = "IPF Trait-Relevant Score (TRS) Across Cell Types"
  )

ggsave(plot = p_trs, filename = paste0(outDir, "/trs_score_violin.png"), width = 10, height = 7, dpi = 300)

# Violin plot of TRS Score by Cell Subtype, faceted by IPF
p_trs_facet <- ggplot(pagwas@meta.data, aes(x = Cell_Subtype, y = scPagwas.downTRS.Score2, color = Cell_Subtype)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  facet_wrap(~IPF) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t=0.1, b=0.1, r=0.1, l=2, "cm"),
    legend.position = "none"
  ) +
  labs(
    x = "Cell Type",
    y = "TRS Score",
    title = "IPF Trait-Relevant Score (TRS) Across Cell Types (Faceted by IPF)"
  )

ggsave(plot = p_trs_facet, filename = paste0(outDir, "/trs_score_violin_facet.png"), width = 10, height = 7, dpi = 300)

# Plot Violin plots of the top 10 genes with the highest PCC between IPF and non-IPF
top_genes <- cor_df$genes[1:10]


for (gene in top_genes) {
  png(paste0(outDir, "/PCC_", gene, "_violin.png"), width = 10, height = 7, units = "in", res = 300)
  print(
    VlnPlot(object = pagwas, features = gene, split.by = 'IPF') +
      theme(
        plot.margin = margin(t = 0.5, b = 0.5, r = 0.5, l = 1, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  )
  dev.off()
}


# Pathway Gene expression plots ---
library(KEGGREST)
pathways.geneList <- keggLink("pathway", "hsa")
pathways.pathwayList <- keggList("pathway", "hsa")

path_df = data.frame(
    gene_id    = gsub("hsa:", "", names(pathways.list)),
    pathway_id = gsub("path:", "", pathways.list)
)

mdm5_kegg <- gg_dot |> filter(celltypes == "mdm_5") |> filter(logrankvalue > 0)

pathway_genes <- path_df |> filter(pathway_id %in% "hsa04350") |> pull(gene_id) |> unique()

library(org.Hs.eg.db)

# Convert Entrez IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = pathway_genes, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

# Filter out any NA values
gene_symbols <- gene_symbols[!is.na(gene_symbols)]