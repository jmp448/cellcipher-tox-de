import scanpy as sc
import pandas as pd
import numpy as np

scvi_embedded_adata_loc = snakemake.input[0]
umap_embedding_loc = snakemake.output[0]
umap_embedded_adata_loc = snakemake.output[1]
cluster_assignment_loc = snakemake.output[2]

# Load data
adata = sc.read_h5ad(scvi_embedded_adata_loc)

# Get neighbors graph
sc.pp.neighbors(adata, use_rep="X_scVI")

# Get umap embedding
sc.tl.umap(adata)
umap_embedding = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs.index)
umap_embedding.to_csv(umap_embedding_loc, sep="\t")

# Get clusters
sc.tl.leiden(adata, resolution=0.25)
cluster_assignments = adata.obs[['leiden']]
cluster_assignments.to_csv(cluster_assignment_loc, sep="\t")

# Save full object
adata.write_h5ad(umap_embedded_adata_loc)
