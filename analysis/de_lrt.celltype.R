library(DESeq2)
library(tidyverse)

# ECTODERM
# Assess overdispersion estimation
ectoderm_dds <- readRDS("results/tox_de.lrt.ectoderm.rds")
png("figs/ectoderm_dispersion.png")
plotDispEsts(ectoderm_dds)
dev.off()

# Plot log fold-changes
ectoderm_result <- results(ectoderm_dds, contrast=c("toxicity", "yes", "no"))
png("figs/ectoderm_wald.png")
plotMA(ectoderm_result, ylim=c(-2, 2))
dev.off()

# Apply celltype-specific shrinkage and re-visualize log fold-changes 
ectoderm_result_shrunk <- lfcShrink(ectoderm_dds, contrast=c("toxicity", "yes", "no"), type="normal")
png("figs/ectoderm_wald_shrunk.png")
plotMA(ectoderm_result_shrunk, ylim=c(-2, 2))
dev.off()

# MESODERM
# Assess overdispersion estimation
mesoderm_dds <- readRDS("results/tox_de.lrt.mesoderm.rds")
png("figs/mesoderm_dispersion.png")
plotDispEsts(mesoderm_dds)
dev.off()

# Plot log fold-changes
mesoderm_result <- results(mesoderm_dds, contrast=c("toxicity", "yes", "no"))
png("figs/mesoderm_wald.png")
plotMA(mesoderm_result, ylim=c(-2, 2))
dev.off()

# Apply celltype-specific shrinkage and re-visualize log fold-changes 
mesoderm_result_shrunk <- lfcShrink(mesoderm_dds, contrast=c("toxicity", "yes", "no"), type="normal")
png("figs/mesoderm_wald_shrunk.png")
plotMA(mesoderm_result_shrunk, ylim=c(-2, 2))
dev.off()

# ENDODERM
# Assess overdispersion estimation
endoderm_dds <- readRDS("results/tox_de.lrt.endoderm.rds")
png("figs/endoderm_dispersion.png")
plotDispEsts(endoderm_dds)
dev.off()

# Plot log fold-changes
endoderm_result <- results(endoderm_dds, contrast=c("toxicity", "yes", "no"))
png("figs/endoderm_wald.png")
plotMA(endoderm_result, ylim=c(-2, 2))
dev.off()

# Apply celltype-specific shrinkage and re-visualize log fold-changes 
endoderm_result_shrunk <- lfcShrink(endoderm_dds, contrast=c("toxicity", "yes", "no"), type="normal")
png("figs/endoderm_wald_shrunk.png")
plotMA(endoderm_result_shrunk, ylim=c(-2, 2))
dev.off()

# PLURIPOTENT
# Assess overdispersion estimation
pluripotent_dds <- readRDS("results/tox_de.lrt.pluripotent.rds")
png("figs/pluripotent_dispersion.png")
plotDispEsts(pluripotent_dds)
dev.off()

# Plot log fold-changes
pluripotent_result <- results(pluripotent_dds, contrast=c("toxicity", "yes", "no"))
png("figs/pluripotent_wald.png")
plotMA(pluripotent_result, ylim=c(-2, 2))
dev.off()

# Apply celltype-specific shrinkage and re-visualize log fold-changes 
pluripotent_result_shrunk <- lfcShrink(pluripotent_dds, contrast=c("toxicity", "yes", "no"), type="normal")
png("figs/pluripotent_wald_shrunk.png")
plotMA(pluripotent_result_shrunk, ylim=c(-2, 2))
dev.off()

# VOLCANO PLOT
ectoderm_result_tbl <- as_tibble(ectoderm_result, rownames="gene") %>%
mutate(celltype="ectoderm")
mesoderm_result_tbl <- as_tibble(mesoderm_result, rownames="gene") %>%
mutate(celltype="mesoderm")
endoderm_result_tbl <- as_tibble(endoderm_result, rownames="gene") %>%
mutate(celltype="endoderm")
pluripotent_result_tbl <- as_tibble(pluripotent_result, rownames="gene") %>%
mutate(celltype="pluripotent")

merged_result_tbl <- bind_rows(ectoderm_result_tbl, mesoderm_result_tbl, endoderm_result_tbl, pluripotent_result_tbl) %>%
mutate(cols=if_else((padj <= 0.05) & log2FoldChange > 0, "red", "gray")) %>%
mutate(cols=if_else((padj <= 0.05) & log2FoldChange < 0, "blue", cols))

png("figs/volcano.png", width=1200, height=400)
ggplot(merged_result_tbl, aes(x=log2FoldChange, y=-log10(pvalue), color=cols)) +
geom_point() +
theme_classic() +
scale_color_identity() +
facet_grid(cols=vars(celltype))
dev.off()

# MASH PREP
betas <- merged_result_tbl %>%
dplyr::select(gene,log2FoldChange,celltype) %>%
pivot_wider(id_cols=gene, names_from=celltype, values_from=log2FoldChange)

ses <- merged_result_tbl %>%
dplyr::select(gene,lfcSE,celltype) %>%
pivot_wider(id_cols=gene, names_from=celltype, values_from=lfcSE)

save(betas, ses, file="results/tox_de.mash_inputs.RData")