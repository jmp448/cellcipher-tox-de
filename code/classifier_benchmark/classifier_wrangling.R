library(tidyverse)

endoderm_de <- readRDS("results/endoderm_de.RData")
mesoderm_de <- readRDS("results/mesoderm_de.RData")
ectoderm_de <- readRDS("results/ectoderm_de.RData")
pluri_de <- readRDS("results/pluripotent_de.RData")

treatments <- levels(colData(endoderm_de)$treatment)[-c(1)] # list treatments besides DMSO

# Create a tibble of all coefficients and p-values for each treatment
pull_coefs <- function(t, deseq) {
  as_tibble(results(deseq, contrast=c("treatment", t, "DMSO")), rownames="gene") %>%
    mutate(treatment=t)
}

endoderm_all <- bind_rows(map_df(treatments, pull_coefs, deseq=endoderm_de))
mesoderm_all <- bind_rows(map_df(treatments, pull_coefs, deseq=mesoderm_de))
ectoderm_all <- bind_rows(map_df(treatments, pull_coefs, deseq=ectoderm_de))
pluri_all <- bind_rows(map_df(treatments, pull_coefs, deseq=pluri_de))

all_treatment_effects <- bind_rows(mutate(endoderm_all, celltype="endoderm"), 
                                   mutate(mesoderm_all, celltype="mesoderm"),
                                   mutate(ectoderm_all, celltype="ectoderm"),
                                   mutate(pluri_all, celltype="pluripotent"))
write_tsv(all_treatment_effects, "results/all_treatment_effects.tsv")

# feature selection
### in each celltype, get 500 genes with the most highly variable effects across treatments
highly_variable_effects <- all_treatment_effects %>%
  group_by(celltype, gene) %>%
  summarize(var_lfc=var(log2FoldChange)) %>%
  group_by(celltype) %>%
  arrange(desc(var_lfc)) %>%
  slice_head(n=500) %>%
  dplyr::select(celltype, gene)

# load the metadata
toxicity_labels <- read_tsv("data/combined_metadata_merged_celltypes.tsv") %>%
  dplyr::select(treatment, toxicity) %>%
  mutate(toxic=if_else(toxicity=="yes", T, F)) %>%
  dplyr::select(treatment, toxic) %>%
  distinct()

classifier_tidied <- all_treatment_effects %>%
  dplyr::select(celltype, gene, treatment, stat) %>%
  inner_join(highly_variable_effects, by=c("celltype", "gene")) %>%
  unite(ctg, celltype, gene, sep="_") %>%
  pivot_wider(id_cols=treatment, names_from=ctg, values_from=stat) %>%
  left_join(toxicity_labels, by=c("treatment")) %>%
  mutate(toxic=factor(toxic)) %>%
  write_tsv("results/tidy_classifier_input.tsv")

