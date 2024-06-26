import scanpy as sc
import anndata as ad

adatas=snakemake.input
full_h5ad = snakemake.output['full']

# Create anndata object
adata_1 = sc.read_h5ad(adatas[0])
adata_2 = sc.read_h5ad(adatas[1])
adata_3 = sc.read_h5ad(adatas[2])
adata_4 = sc.read_h5ad(adatas[3])
adata_5 = sc.read_h5ad(adatas[4])
adata_6 = sc.read_h5ad(adatas[5])
adata_7 = sc.read_h5ad(adatas[6])

adata_full = ad.concat([adata_1, adata_2, adata_3, adata_4, adata_5, adata_6, adata_7])

adata_full.write_h5ad(filename=full_h5ad)
