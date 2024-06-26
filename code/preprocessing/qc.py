import scanpy as sc
import numpy as np
import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix
import seaborn as sns
import matplotlib.pyplot as plt

full_h5ad_loc = snakemake.input['full_adata']
filtered_h5ad_loc = snakemake.output['filtered_adata']
metadata_loc = snakemake.output['metadata']


# Create anndata object
adata = sc.read_h5ad(full_h5ad_loc)

# Filter anndata object
adata = adata[adata.obs['sample_id'] != "Undetermined"]
adata = adata[~adata.obs['donor_id'].isin(['unassigned', 'doublet'])]
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.pct_counts_mt < 15, :]
adata = adata[adata.obs.total_counts < 40000, :]

# Add treatment labels
# treatment_map = {1:'DMSO', 2:'DMSO', 3:'Difloxacin', 4:'Difloxacin', 5:'Ambrisentan', 6:'Ambrisentan', 7:'BIA', 8:'BIA', 9:'Fialuridine', 10:'Fialuridine', 11:'Benfluorex_DMSO', 12:'Benfluorex_DMSO', 13:'Benfluorex_H2O', 14:'Benfluorex_H2O', 15:'Pergolide', 16:'Pergolide'}
adata.obs['treatment'] = [sample.split('-')[-1] for sample in adata.obs['sample_id']]

adata.write_h5ad(filtered_h5ad_loc)
adata.obs.to_csv(metadata_loc, sep='\t')
