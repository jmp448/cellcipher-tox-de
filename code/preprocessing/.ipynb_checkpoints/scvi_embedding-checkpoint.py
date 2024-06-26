import torch
import scanpy as sc
import pandas as pd
import numpy as np
import scvi
import random

random.seed(1234)

adata_loc = snakemake.input[0]
model_loc = snakemake.params[0]
embedding_loc = snakemake.output[0]
adata_embedded_loc = snakemake.output[1]

print(torch.cuda.is_available())
print(torch.cuda.device_count())

# Load and tidy data
adata = sc.read_h5ad(adata_loc)
adata.obs = adata.obs.astype({"donor_id":'category', "sample_id":'category'})

# Set up anndata object for scVI
scvi.model.SCVI.setup_anndata(adata, layer="counts",
                              categorical_covariate_keys=['sample_id'])

# Train linearly decoded VAE model
mod = scvi.model.SCVI(adata, n_latent=50, gene_likelihood = 'nb')
mod.train(max_epochs=100, plan_kwargs={'lr':5e-3}, check_val_every_n_epoch=10)
mod.save(model_loc)

# Save results
adata.obsm["X_scVI"] = mod.get_latent_representation()
scvi_embedding = pd.DataFrame(adata.obsm['X_scVI'], index=adata.obs.index)
scvi_embedding.to_csv(embedding_loc, sep="\t")
adata.write_h5ad(adata_embedded_loc)
