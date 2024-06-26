library(DESeq2)
library(edgeR)
library(tidyverse)
library(Matrix)
library(Matrix.utils)
library(matrixStats)
library(vroom)
library(ggrepel)
library(RColorBrewer)
library(viridis)

load("data/expression_filtered.RData")
sample_metadata_cv <- filter(sample_metadata_filt, cv_tox != "n/a")
expression.mat.filtered_cv <- expression.mat.filtered_filt[,sample_metadata_cv$pseudobulk_sample]

ct_de_cv_toxicity <- function(t) {
  dds_metadata_cv <- sample_metadata_cv %>%
    filter(celltype == !!t) %>%
    column_to_rownames("pseudobulk_sample")
  dds_metadata_cv$donor_id <- factor(dds_metadata_cv$donor_id)
  dds_metadata_cv$cv_tox <- factor(dds_metadata_cv$cv_tox, levels=c("yes", "no"))
  dds_counts_cv <- expression.mat.filtered_cv[,rownames(dds_metadata_cv)]

  dds_cv <- DESeqDataSetFromMatrix(countData = dds_counts_cv,
                                colData = dds_metadata_cv,
                                design = ~ donor_id + cv_tox)
  dds_cv <- DESeq(dds_cv)

  as_tibble(results(dds_cv, contrast=c("cv_tox", "yes", "no"))@listData) %>%
    mutate(celltype=!!t)
}

celltypes <- unique(sample_metadata_cv$celltype)
tox_comp_out <- bind_rows(map_df(celltypes, ct_de_cv_toxicity))
tox_comp_out$gene <- rep(rownames(expression.mat.filtered_cv), length(celltypes))
tox_comp_out$nlog10p <- -log10(tox_comp_out$pvalue)

tox_comp_out$de <- if_else(tox_comp_out$nlog10p >= 5 & tox_comp_out$log2FoldChange <= -0.25, "DOWN", "Not DE")
tox_comp_out$de <- if_else(tox_comp_out$nlog10p >= 5 & tox_comp_out$log2FoldChange >= 0.25, "UP", tox_comp_out$de)
tox_comp_out$de <- factor(tox_comp_out$de, levels=c("DOWN", "Not DE", "UP"))
tox_comp_out$label <- if_else(tox_comp_out$de != "Not DE", tox_comp_out$gene, NA_character_)

ggplot(tox_comp_out, aes(x=log2FoldChange, y=nlog10p, color=de, label=label)) +
  geom_point() +
  facet_grid(cols=vars(celltype)) +
  theme_classic(base_size=12) +
  scale_color_manual(values=c("#B51700", "#BBBBBB", "#00A1FF")) +
  geom_text_repel() +
  ggtitle("cv Toxicity DE")
write_tsv(tox_comp_out, "temp/de_comparison_cv_toxicity.tsv")

all_de_genes <- tox_comp_out %>%
  filter(padj <= 0.05) %>%
  select(celltype, gene, padj, log2FoldChange) %>%
  distinct() %>%
  arrange(padj) %>%
  write_tsv("temp/de_genes_batches_12456_cv.tsv")

celltype_upregulated_genes <- tox_comp_out %>%
  filter(log2FoldChange > 0) %>%
  filter(padj <= 0.1) %>%
  select(celltype, gene)

## GSEA
gmt_file <- "/project2/gilad/jpopp/ebQTL/data/gene_sets/c5.go.bp.v2022.1.Hs.symbols.gmt"
gmt_lines <- readLines(gmt_file)
gmt_list <- lapply(gmt_lines, function(x) unlist(strsplit(x, "\t")))
gmt_list2 <- lapply(gmt_list, function(x) list(x[1], x[2], x[3:length(x)]))
gmt_df <- as_tibble(do.call("rbind", gmt_list2)) %>%
  mutate(across(V1:V2, unlist)) %>%
  dplyr::rename(geneset=V1, link=V2, genes=V3) %>%
  filter(str_sub(geneset, 1, 4) == "GOBP") # filter to BP 
gmt_mat <- gmt_df %>%
  dplyr::select(-c(link)) %>%
  unnest(genes) %>%
  filter(genes %in% unique(full_out$gene)) %>%
  pivot_wider(names_from = genes, values_from = genes, values_fn = length, values_fill = 0) %>%
  column_to_rownames("geneset")
valid_genesets <- (10 <= rowSums(gmt_mat)) & (rowSums(gmt_mat) < 200)
gmt_mat <- t(as.matrix(gmt_mat[which(valid_genesets),]))

## cv toxicity
toxicity_map <- vroom("data/tox_toxicity.txt", col_names=c("treatment", "toxicity")) %>%
  distinct()

treatment_dict <- treatment_dict %>%
  select(-c(treatment)) %>%
  dplyr::rename(treatment = treatment_name) %>%
  left_join(toxicity_map, by="treatment")


