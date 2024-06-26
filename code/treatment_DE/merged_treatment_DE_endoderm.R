library(DESeq2)
library(tidyverse)
library(BiocParallel)
register(MulticoreParam(40))
print(parallel::detectCores())

load("data/expression_filtered_merged_celltypes.RData")

dds_metadata_filt <- sample_metadata_filt %>%
  filter(celltype == "endoderm") %>%
  column_to_rownames("pseudobulk_sample") 
dds_metadata_filt$donor_id <- factor(dds_metadata_filt$donor_id)
dds_metadata_filt$treatment <- factor(dds_metadata_filt$treatment)
dds_metadata_filt$treatment <- relevel(dds_metadata_filt$treatment, ref = "DMSO")
dds_counts_filt <- expression.mat.filtered_filt[,rownames(dds_metadata_filt)]

dds_filt <- DESeqDataSetFromMatrix(countData = dds_counts_filt,
                                   colData = dds_metadata_filt,
                                   design = ~ donor_id + treatment)
dds_filt <- DESeq(dds_filt, parallel=T)

saveRDS(dds_filt, "results/endoderm_de.RData")
