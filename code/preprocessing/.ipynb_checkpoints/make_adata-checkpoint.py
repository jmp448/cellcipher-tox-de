import scanpy as sc
import numpy as np
import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix

aggregator=snakemake.input['aggregator']
combined_reads_dir = snakemake.params['prefix']
full_h5ad = snakemake.output['outfile']

# Create anndata object
adata = sc.read_10x_mtx(combined_reads_dir)
adata.X = csr_matrix(adata.X)

# Add sample metadata
keys = [ s[s.find("-")+1:] for s in list(adata.obs.index)]
keys = pd.to_numeric(keys)-1
adata.obs['key']=keys

sample_key = pd.read_csv(aggregator)[['sample_id']].reset_index().rename(columns={'index': 'key'})

sample_metadata = adata.obs.reset_index().rename(columns={'index': 'cell'}).merge(sample_key, on='key', how='left').set_index('cell')
adata.obs = sample_metadata.drop(columns='key')

# Add vireo metadata
adata.obs['sample_barcode'] = [b[:b.find("-")] for b in list(adata.obs.index)]

vireo_dfs = []
for s in list(np.unique(adata.obs['sample_id'])):
    v = pd.read_csv(f"data_071023/vireo/{s}/donor_ids.tsv", sep="\t")
    v['sample_id'] = s
    v['sample_barcode'] = [b[:b.find("-")] for b in list(v['cell'])]
    vireo_dfs.append(v)

vireo_df = pd.concat(vireo_dfs).drop(columns='cell')
demultiplex_metadata = adata.obs.reset_index().rename(columns={'index': 'cell'}).merge(vireo_df, on=['sample_barcode', 'sample_id'], how='left').set_index('cell')

adata.obs = demultiplex_metadata

adata.write_h5ad(filename=full_h5ad)
