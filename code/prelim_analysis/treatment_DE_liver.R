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
sample_metadata_liver <- filter(sample_metadata_filt, liver_tox != "n/a")
expression.mat.filtered_liver <- expression.mat.filtered_filt[,sample_metadata_liver$pseudobulk_sample]

ct_de_liver_toxicity <- function(t) {
  dds_metadata_liver <- sample_metadata_liver %>%
    filter(celltype == !!t) %>%
    column_to_rownames("pseudobulk_sample")
  dds_metadata_liver$donor_id <- factor(dds_metadata_liver$donor_id)
  dds_metadata_liver$liver_tox <- factor(dds_metadata_liver$liver_tox, levels=c("yes", "no"))
  dds_counts_liver <- expression.mat.filtered_liver[,rownames(dds_metadata_liver)]

  dds_liver <- DESeqDataSetFromMatrix(countData = dds_counts_liver,
                                colData = dds_metadata_liver,
                                design = ~ donor_id + liver_tox)
  dds_liver <- DESeq(dds_liver)

  as_tibble(results(dds_liver, contrast=c("liver_tox", "yes", "no"))@listData) %>%
    mutate(celltype=!!t)
}

celltypes <- unique(sample_metadata_liver$celltype)
tox_comp_out <- bind_rows(map_df(celltypes, ct_de_liver_toxicity))
tox_comp_out$gene <- rep(rownames(expression.mat.filtered_liver), length(celltypes))
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
  ggtitle("Liver Toxicity DE")
write_tsv(tox_comp_out, "temp/de_comparison_liver_toxicity.tsv")

all_de_genes <- tox_comp_out %>%
  filter(padj <= 0.05) %>%
  select(celltype, gene, padj, log2FoldChange) %>%
  distinct() %>%
  arrange(padj) %>%
  write_tsv("temp/de_genes_batches_12456_liver.tsv")

