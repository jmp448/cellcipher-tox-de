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

toxicity_map <- vroom("data/tox_toxicity.txt") %>%
  distinct()

treatment_dict <- treatment_dict %>%
  select(-c(treatment)) %>%
  dplyr::rename(treatment = treatment_name) %>%
  left_join(toxicity_map, by="treatment")

## PSEUDOBULK DATA
pseudobulk_exp_1 <- vroom("/project/gilad/jpopp/cellcipher/data/Tox1-full/pseudobulk_expression.tsv") %>%
  dplyr::rename(gene=`...1`)
pseudobulk_exp_2 <- vroom("/project/gilad/jpopp/cellcipher/data/Tox2-full/pseudobulk_expression.tsv") %>%
  dplyr::rename(gene=`...1`)
pseudobulk_exp_4 <- vroom("/project/gilad/jpopp/cellcipher/data/Tox4-full/pseudobulk_expression.tsv") %>%
  dplyr::rename(gene=`...1`)
pseudobulk_exp_5 <- vroom("/project/gilad/jpopp/cellcipher/data/Tox5-full/pseudobulk_expression.tsv") %>%
  dplyr::rename(gene=`...1`)
pseudobulk_exp_6 <- vroom("/project/gilad/jpopp/cellcipher/data/Tox6-full/pseudobulk_expression.tsv") %>%
  dplyr::rename(gene=`...1`)

pseudobulk_exp <- inner_join(pseudobulk_exp_1, pseudobulk_exp_2, by="gene") %>%
  inner_join(pseudobulk_exp_4, by="gene") %>%
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

sample_metadata_filt <- filter(sample_metadata, !treatment %in% c("Benfluorex_DMSO", "DMSO"))
pseudobulk_exp_filt <- select(pseudobulk_exp, c(gene, sample_metadata_filt$pseudobulk_sample))

expression.mat_filt <- pseudobulk_exp_filt %>% column_to_rownames("gene") %>% as.matrix
samples_filt <- colnames(expression.mat_filt)
pcgs_filt <- rownames(expression.mat_filt) %in% filtered_gtf$hgnc
gene_vars_filt <- rowVars(expression.mat_filt) > 0
gene_nonzeros_filt <- rowMedians(expression.mat_filt) > 10
gene_keepers_filt <- rownames(expression.mat_filt)[pcgs_filt & gene_vars_filt & gene_nonzeros_filt]
expression.mat.filtered_filt <- expression.mat_filt[gene_keepers_filt,]

save(sample_metadata_filt, pseudobulk_exp_filt, expression.mat.filtered_filt,  
     file="data/expression_filtered.RData")