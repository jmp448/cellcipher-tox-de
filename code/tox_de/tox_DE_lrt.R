library(DESeq2)
library(tidyverse)
library(BiocParallel)

load("data/training/tox_de_training_data.RData")

dds_metadata <- training_metadata %>%
  dplyr::select(pseudobulk_sample, donor_id, batch, celltype, toxicity, treatment) %>%
  column_to_rownames("pseudobulk_sample") 
dds_metadata$donor_id <- factor(dds_metadata$donor_id)
dds_metadata$batch <- factor(dds_metadata$batch)
dds_metadata$celltype <- factor(dds_metadata$celltype)
dds_metadata$toxicity <- factor(dds_metadata$toxicity)
dds_counts <- training_pseudobulk %>%
  column_to_rownames("gene") %>% as.matrix
dds_counts <- dds_counts[,rownames(dds_metadata)]

dds <- DESeqDataSetFromMatrix(countData = dds_counts,
                                   colData = dds_metadata,
                                   design = ~ donor_id + batch + celltype + toxicity + celltype * toxicity)
system.time(
  deseq_full <- DESeq(dds, test="LRT",
                      full= ~ donor_id + batch + celltype + toxicity + celltype * toxicity,
                      reduced = ~ donor_id + batch + celltype,
                      parallel=T)
)
print("finished model fit")
saveRDS(deseq_full, "results/tox_de.lrt.rds")