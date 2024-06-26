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

## TREATMENT METADATA
treatment_dict_1 <- vroom("/project/gilad/jpopp/cellcipher/data/Tox1-full/treatment_dict.tsv")
treatment_dict_2 <- vroom("/project/gilad/jpopp/cellcipher/data/Tox2-full/treatment_dict.tsv")
treatment_dict_4 <- vroom("/project/gilad/jpopp/cellcipher/data/Tox4-full/treatment_dict.tsv") %>%
  mutate(treatment=as.character(treatment))
treatment_dict_5 <- vroom("/project/gilad/jpopp/cellcipher/data/Tox5-full/treatment_dict.tsv")
treatment_dict_6 <- vroom("/project/gilad/jpopp/cellcipher/data/Tox6-full/treatment_dict.tsv") %>%
  mutate(treatment=as.character(treatment))

treatment_dict <- dplyr::bind_rows(treatment_dict_1,
                                   treatment_dict_2,
                                   treatment_dict_4,
                                   treatment_dict_5,
                                   treatment_dict_6)

toxicity_map <- vroom("data/tox_toxicity.txt", col_names=c("treatment", "toxicity")) %>%
  distinct()

treatment_dict <- treatment_dict %>%
  select(-c(treatment)) %>%
  dplyr::rename(treatment = treatment_name) %>%
  left_join(toxicity_map, by="treatment")

## PSEUDOBULK DATA
pseudobulk_exp_1 <- vroom("/project/gilad/jpopp/cellcipher/data/Tox1-full/pseudobulk_expression.tsv") %>%
  dplyr::rename(gene=`...1`)
pseudobulk_exp_4 <- vroom("/project/gilad/jpopp/cellcipher/data/Tox4-full/pseudobulk_expression.tsv") %>%
  dplyr::rename(gene=`...1`)
pseudobulk_exp_5 <- vroom("/project/gilad/jpopp/cellcipher/data/Tox5-full/pseudobulk_expression.tsv") %>%
  dplyr::rename(gene=`...1`)
pseudobulk_exp_6 <- vroom("/project/gilad/jpopp/cellcipher/data/Tox6-full/pseudobulk_expression.tsv") %>%
  dplyr::rename(gene=`...1`)

pseudobulk_exp <- inner_join(pseudobulk_exp_1, pseudobulk_exp_4, by="gene") %>%
  inner_join(pseudobulk_exp_5, by="gene") %>%
  inner_join(pseudobulk_exp_6, by="gene")

sample_metadata <- tibble("pseudobulk_sample"=colnames(pseudobulk_exp)[2:ncol(pseudobulk_exp)]) %>%
  separate(pseudobulk_sample, into=c("sample_id", "donor_id", "celltype"), sep="_", remove=F) %>%
  left_join(treatment_dict, by="sample_id")

write_tsv(pseudobulk_exp, "data/combined_pseudobulk.tsv")
write_tsv(sample_metadata, "data/combined_metadata.tsv")

pseudobulk_exp <- vroom("data/combined_pseudobulk.tsv")
sample_metadata <- vroom("data/combined_metadata.tsv")

# Filter & Normalize
filtered_gtf <- vroom("/project2/gilad/jpopp/ebQTL/data/gencode/gencode.hg38.filtered.gtf")

invnorm_transform <- function(x) {
  qqnorm(x, plot.it=F)$x
}

# Take the sum over replicates (? skipping this ?)
# pseudobulk_exp <- pseudobulk_exp %>%
#   pivot_longer(!gene, names_to="pseudobulk_sample", values_to="counts") %>%
#   left_join(sample_metadata, by="pseudobulk_sample") %>%
#   unite(aggregate_sample, treatment, donor_id, celltype, sep="_") %>%
#   group_by(aggregate_sample, gene) %>%
#   summarise(aggregate_counts=sum(counts)) %>%
#   pivot_wider(id_cols=gene, names_from=aggregate_sample, values_from=aggregate_counts)
# sample_metadata <- sample_metadata %>%
#   unite(pseudobulk_sample, treatment, donor_id, celltype, sep="_", remove=F) %>%
#   select(-sample_id) %>%
#   distinct()

# Filter to protein coding genes with nonzero variance, nonzero median
expression.mat <- pseudobulk_exp %>% column_to_rownames("gene") %>% as.matrix
samples <- colnames(expression.mat)
pcgs <- rownames(expression.mat) %in% filtered_gtf$hgnc
gene_vars <- rowVars(expression.mat) > 0
gene_nonzeros <- rowMedians(expression.mat) > 10
gene_keepers <- rownames(expression.mat)[pcgs & gene_vars & gene_nonzeros]
expression.mat.filtered <- expression.mat[gene_keepers,]

# Run TMM normalization
# 1. TMM normalization to get log TMM-normalized counts
expression.mat.tmm <- DGEList(counts=expression.mat.filtered) %>% calcNormFactors(method="TMM")
expression.mat.logcpm <- cpm(expression.mat.tmm, log=T)
# 2. normalize each gene to have zero mean, unit variance across individuals
expression.mat.std <- apply(expression.mat.logcpm, 1, invnorm_transform) %>% t
# 3. re-insert gene names and convert to tibble
rownames(expression.mat.std) <- gene_keepers
colnames(expression.mat.std) <- samples

# Check for any outliers w PCA
expression.pca <- princomp(expression.mat.std)

loadings <- as_tibble(expression.pca$loadings[,1:10], rownames="pseudobulk_sample") %>%
  left_join(sample_metadata, by="pseudobulk_sample") %>%
  rename(PC1=Comp.1, PC2=Comp.2, PC3=Comp.3, PC4=Comp.4, PC5=Comp.5,
         PC6=Comp.6, PC7=Comp.7, PC8=Comp.8, PC9=Comp.9, PC10=Comp.10)

ggplot(loadings, aes(x=PC1, y=PC2, color=celltype)) +
  geom_point() +
  theme_classic(base_size=15)
# ggplot(loadings, aes(x=PC1, y=PC5, color=treatment)) +
#   geom_point() +
#   # scale_color_manual(values=color_palette) + 
#   theme_classic(base_size=15)
# ggplot(loadings, aes(x=treatment, y=PC5, fill=treatment)) +
#   geom_violin() + 
#   theme_classic(base_size=15) +
#   # scale_fill_manual(values=color_palette) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   
# adj.r2 <- data.frame(matrix(NA, nrow=3, ncol=10)) 
# rownames(adj.r2) <- c("donor_id", "treatment", "celltype")
# colnames(adj.r2) <- paste0("PC", seq(10))
# for (var in rownames(adj.r2)) {
#   for (pc in colnames(adj.r2)) {
#     adj.r2[var, pc] <- summary(lm(paste0(pc, " ~ ", var), loadings))$adj.r.squared
#   }
# }
# as_tibble(adj.r2, rownames="variable") %>%
#   pivot_longer(!variable, names_to="component", values_to="R2") %>%
#   mutate(component=factor(component, levels=paste0("PC", seq(10)))) %>%
#   ggplot(aes(x=variable, y=component, fill=R2)) +
#   geom_tile() +
#   scale_fill_viridis() + 
#   theme_classic(base_size=15)
# 
# # DE genes wrt PC5
# pc5_scores <- expression.pca$scores[,"Comp.5"]
# hist(pc5_scores)
# pc5_up <- as_tibble(head(pc5_scores, n=100), rownames="gene") %>% dplyr::arrange(desc(value))
# pc5_down <- as_tibble(tail(pc5_scores, n=100), rownames="gene") %>% dplyr::arrange(value)
# write_tsv(pc5_up, "temp/Tox1_pc5_upregulated.tsv")
# write_tsv(pc5_down, "temp/Tox1_pc5_downregulated.tsv")
# 
# # How much variance is explained 
# # Perform a differential expression analysis 
# t <- "mesoderm"
# treatments <- c("Difloxacin", "Ambrisentan", "BIA",
#                 "Fialuridine","Benfluorex_DMSO","Benfluorex_H2O","Pergolide")

# ct_de <- function(t) {
#   dds_metadata <- sample_metadata %>%
#     filter(celltype == !!t) %>%
#     column_to_rownames("pseudobulk_sample")
#   dds_metadata$donor_id <- factor(dds_metadata$donor_id)
#   dds_metadata$treatment <- factor(dds_metadata$treatment, levels=c("DMSO", "Difloxacin", "Ambrisentan", "BIA",
#                                                                     "Fialuridine","Benfluorex_DMSO","Benfluorex_H2O","Pergolide"))
#   dds_counts <- expression.mat.filtered[,rownames(dds_metadata)]
#   
#   dds <- DESeqDataSetFromMatrix(countData = dds_counts,
#                                 colData = dds_metadata,
#                                 design = ~ donor_id + treatment)
#   dds <- DESeq(dds)
#   
#   difloxacin <- as_tibble(results(dds, contrast=c("treatment", "Difloxacin", "DMSO"))@listData) %>% mutate(treatment="Difloxacin")
#   ambrisentan <- as_tibble(results(dds, contrast=c("treatment", "Ambrisentan", "DMSO"))@listData) %>% mutate(treatment="Ambrisentan")
#   bia <- as_tibble(results(dds, contrast=c("treatment", "BIA", "DMSO"))@listData) %>% mutate(treatment="BIA")
#   fialuridine <- as_tibble(results(dds, contrast=c("treatment", "Fialuridine", "DMSO"))@listData) %>% mutate(treatment="Fialuridine")
#   benfluorex_dmso <- as_tibble(results(dds, contrast=c("treatment", "Benfluorex_DMSO", "DMSO"))@listData) %>% mutate(treatment="Benfluorex_DMSO")
#   benfluorex_h2o <- as_tibble(results(dds, contrast=c("treatment", "Benfluorex_H2O", "DMSO"))@listData) %>% mutate(treatment="Benfluorex_H2O")
#   pergolide <- as_tibble(results(dds, contrast=c("treatment", "Pergolide", "DMSO"))@listData) %>% mutate(treatment="Pergolide")
#   
#   ct <- bind_rows(difloxacin, ambrisentan, bia, fialuridine, benfluorex_dmso, benfluorex_h2o, pergolide) %>% mutate(celltype=!!t)
# }
# 
# celltypes <- unique(sample_metadata$celltype)
# full_out <- bind_rows(map_df(celltypes, ct_de))
# full_out$gene <- rep(rownames(expression.mat.filtered), length(celltypes)*7)
# full_out$nlog10p <- -log10(full_out$pvalue)
# 
# write_tsv(full_out, "temp/tox1_treatment_DE.tsv")
# 
# full_out <- read_tsv("temp/tox1_treatment_DE.tsv")
# ggplot(full_out, aes(x=log2FoldChange, y=nlog10p)) +
#   geom_point() +
#   facet_grid(rows=vars(treatment), cols=vars(celltype)) +
#   theme_classic(base_size=15)
# 
# # Ambrisentan
# t1 <- filter(full_out, treatment == "Ambrisentan")
# t1$de <- if_else(t1$nlog10p >= 10 & t1$log2FoldChange <= -0.5, "DOWN", "Not DE")
# t1$de <- if_else(t1$nlog10p >= 10 & t1$log2FoldChange >= 0.5, "UP", t1$de)
# t1$de <- factor(t1$de, levels=c("DOWN", "Not DE", "UP"))
# t1$label <- if_else(t1$de != "Not DE", t1$gene, NA_character_)
# 
# ggplot(t1, aes(x=log2FoldChange, y=nlog10p, color=de, label=label)) +
#   geom_point() +
#   facet_grid(cols=vars(celltype)) +
#   theme_classic(base_size=15) +
#   scale_color_manual(values=c("#B51700", "#BBBBBB", "#00A1FF")) +
#   geom_text_repel() +
#   ggtitle("Ambrisentan DE")
# 
# # BIA
# t2 <- filter(full_out, treatment == "BIA")
# t2$de <- if_else(t2$nlog10p >= 10 & t2$log2FoldChange <= -0.5, "DOWN", "Not DE")
# t2$de <- if_else(t2$nlog10p >= 10 & t2$log2FoldChange >= 0.5, "UP", t2$de)
# t2$de <- factor(t2$de, levels=c("DOWN", "Not DE", "UP"))
# t2$label <- if_else(t2$de != "Not DE", t2$gene, NA_character_)
# 
# ggplot(t2, aes(x=log2FoldChange, y=nlog10p, color=de, label=label)) +
#   geom_point() +
#   facet_grid(cols=vars(celltype)) +
#   theme_classic(base_size=15) +
#   scale_color_manual(values=c("#B51700", "#BBBBBB", "#00A1FF")) +
#   geom_text_repel() +
#   ggtitle("BIA DE")
# 
# # Benfluorex DMSO
# t3 <- filter(full_out, treatment == "Benfluorex_DMSO")
# t3$de <- if_else(t3$nlog10p >= 10 & t3$log2FoldChange <= -1.5, "DOWN", "Not DE")
# t3$de <- if_else(t3$nlog10p >= 10 & t3$log2FoldChange >= 1, "UP", t3$de)
# t3$de <- factor(t3$de, levels=c("DOWN", "Not DE", "UP"))
# t3$label <- if_else(t3$de != "Not DE", t3$gene, NA_character_)
# 
# ggplot(t3, aes(x=log2FoldChange, y=nlog10p, color=de, label=label)) +
#   geom_point() +
#   facet_grid(cols=vars(celltype)) +
#   theme_classic(base_size=15) +
#   scale_color_manual(values=c("#B51700", "#BBBBBB", "#00A1FF")) +
#   geom_text_repel() +
#   ggtitle("Benfluorex DMSO DE")
# 
# # Benfluorex H2O
# t4 <- filter(full_out, treatment == "Benfluorex_H2O")
# t4$de <- if_else(t4$nlog10p >= 10 & t4$log2FoldChange <= -1.5, "DOWN", "Not DE")
# t4$de <- if_else(t4$nlog10p >= 10 & t4$log2FoldChange >= 1.5, "UP", t4$de)
# t4$de <- factor(t4$de, levels=c("DOWN", "Not DE", "UP"))
# t4$label <- if_else(t4$de != "Not DE", t4$gene, NA_character_)
# 
# ggplot(t4, aes(x=log2FoldChange, y=nlog10p, color=de, label=label)) +
#   geom_point() +
#   facet_grid(cols=vars(celltype)) +
#   theme_classic(base_size=15) +
#   scale_color_manual(values=c("#B51700", "#BBBBBB", "#00A1FF")) +
#   geom_text_repel() +
#   ggtitle("Benfluorex H2O DE")
# 
# # Difloxacin
# t5 <- filter(full_out, treatment == "Difloxacin")
# t5$de <- if_else(t5$nlog10p >= 10 & t5$log2FoldChange <= -0.5, "DOWN", "Not DE")
# t5$de <- if_else(t5$nlog10p >= 10 & t5$log2FoldChange >= 0.5, "UP", t5$de)
# t5$de <- factor(t5$de, levels=c("DOWN", "Not DE", "UP"))
# t5$label <- if_else(t5$de != "Not DE", t5$gene, NA_character_)
# 
# ggplot(t5, aes(x=log2FoldChange, y=nlog10p, color=de, label=label)) +
#   geom_point() +
#   facet_grid(cols=vars(celltype)) +
#   theme_classic(base_size=15) +
#   scale_color_manual(values=c("#B51700", "#BBBBBB", "#00A1FF")) +
#   geom_text_repel() +
#   ggtitle("Difloxacin DE")
# 
# # Fialuridine
# t6 <- filter(full_out, treatment == "Fialuridine")
# t6$de <- if_else(t6$nlog10p >= 10 & t6$log2FoldChange <= -1, "DOWN", "Not DE")
# t6$de <- if_else(t6$nlog10p >= 10 & t6$log2FoldChange >= 0.5, "UP", t6$de)
# t6$de <- factor(t6$de, levels=c("DOWN", "Not DE", "UP"))
# t6$label <- if_else(t6$de != "Not DE", t6$gene, NA_character_)
# 
# ggplot(t6, aes(x=log2FoldChange, y=nlog10p, color=de, label=label)) +
#   geom_point() +
#   facet_grid(cols=vars(celltype)) +
#   theme_classic(base_size=15) +
#   scale_color_manual(values=c("#B51700", "#BBBBBB", "#00A1FF")) +
#   geom_text_repel() +
#   ggtitle("Fialuridine DE")
# 
# # Pergolide
# t7 <- filter(full_out, treatment == "Pergolide")
# t7$de <- if_else(t7$nlog10p >= 6 & t7$log2FoldChange <= -0.25, "DOWN", "Not DE")
# t7$de <- if_else(t7$nlog10p >= 6 & t7$log2FoldChange >= 0.25, "UP", t7$de)
# t7$de <- factor(t7$de, levels=c("DOWN", "Not DE", "UP"))
# t7$label <- if_else(t7$de != "Not DE", t7$gene, NA_character_)
# 
# ggplot(t7, aes(x=log2FoldChange, y=nlog10p, color=de, label=label)) +
#   geom_point() +
#   facet_grid(cols=vars(celltype)) +
#   theme_classic(base_size=15) +
#   scale_color_manual(values=c("#B51700", "#BBBBBB", "#00A1FF")) +
#   geom_text_repel() +
#   ggtitle("Pergolide DE")
# 
# DE between non-toxic and toxic 
# sample_metadata_noben <- filter(sample_metadata, treatment != "Benfluorex_DMSO") %>%
#   mutate(toxicity = if_else(treatment %in% c("DMSO", "Ambrisentan"), "non-toxic", "toxic"))
sample_metadata_noben <- filter(sample_metadata, treatment != "Benfluorex_DMSO")
pseudobulk_exp_noben <- select(pseudobulk_exp, c(gene, sample_metadata_noben$pseudobulk_sample))

expression.mat_noben <- pseudobulk_exp_noben %>% column_to_rownames("gene") %>% as.matrix
samples_noben <- colnames(expression.mat_noben)
pcgs_noben <- rownames(expression.mat_noben) %in% filtered_gtf$hgnc
gene_vars_noben <- rowVars(expression.mat_noben) > 0
gene_nonzeros_noben <- rowMedians(expression.mat_noben) > 10
gene_keepers_noben <- rownames(expression.mat_noben)[pcgs_noben & gene_vars_noben & gene_nonzeros_noben]
expression.mat.filtered_noben <- expression.mat_noben[gene_keepers_noben,]

ct_de_toxicity <- function(t) {
  dds_metadata_noben <- sample_metadata_noben %>%
    filter(celltype == !!t) %>%
    column_to_rownames("pseudobulk_sample")
  dds_metadata_noben$donor_id <- factor(dds_metadata_noben$donor_id)
  dds_metadata_noben$toxicity <- factor(dds_metadata_noben$toxicity, levels=c("yes", "no"))
  # dds_metadata_noben$toxicitysig <- factor(dds_metadata_noben$toxicitysig, levels=c("toxicity-group-1", "other"))
  dds_counts_noben <- expression.mat.filtered_noben[,rownames(dds_metadata_noben)]

  dds_noben <- DESeqDataSetFromMatrix(countData = dds_counts_noben,
                                colData = dds_metadata_noben,
                                design = ~ donor_id + toxicity)
  # dds_noben <- DESeqDataSetFromMatrix(countData = dds_counts_noben,
  #                                     colData = dds_metadata_noben,
  #                                     design = ~ donor_id + toxicitysig)
  dds_noben <- DESeq(dds_noben)

  as_tibble(results(dds_noben, contrast=c("toxicity", "yes", "no"))@listData) %>%
    mutate(celltype=!!t)
  # as_tibble(results(dds_noben, contrast=c("toxicitysig", "toxicity-group-1", "other"))@listData) %>%
  #   mutate(celltype=!!t)
}

celltypes <- unique(sample_metadata_noben$celltype)
tox_comp_out <- bind_rows(map_df(celltypes, ct_de_toxicity))
tox_comp_out$gene <- rep(rownames(expression.mat.filtered_noben), length(celltypes))
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
  #ggtitle("Toxicity DE")
# write_tsv(tox_comp_out, "temp/de_comparison_toxicity.tsv")
write_tsv(tox_comp_out, "temp/de_comparison_toxicity.tsv")

all_de_genes <- tox_comp_out %>%
  filter(padj <= 0.05) %>%
  select(celltype, gene, padj, log2FoldChange) %>%
  distinct() %>%
  arrange(padj) %>%
  write_tsv("temp/de_genes_batches_1456.tsv")

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

## Liver toxicity
toxicity_map <- vroom("data/tox_toxicity.txt", col_names=c("treatment", "toxicity")) %>%
  distinct()

treatment_dict <- treatment_dict %>%
  select(-c(treatment)) %>%
  dplyr::rename(treatment = treatment_name) %>%
  left_join(toxicity_map, by="treatment")


