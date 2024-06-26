import glob
import os

# samples = set(glob_wildcards("data_071023/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
# determined_samples = samples.copy()
# determined_samples.remove('Undetermined')

def list_donor_ids(wildcards):
  if wildcards.group == "Tox1-full":
    samples = os.listdir("data/Tox1-full/aligned_reads/")
    ids = ['data/Tox1-full/vireo/' + s + '/donor_ids.tsv' for s in samples]
  elif wildcards.group == "Tox2-full":
    samples = os.listdir("data/Tox2-full/aligned_reads/")
    ids = ['data/Tox2-full/vireo/' + s + '/donor_ids.tsv' for s in samples]
  elif wildcards.group == "Tox3-full":
    samples = [s for s in os.listdir("data/Tox3-full/aligned_reads/") if not s.endswith(".mro")]
    ids = ['data/Tox3-full/vireo/' + s + '/donor_ids.tsv' for s in samples]
  elif wildcards.group == "Tox4-full":
    samples = os.listdir("data/Tox4-full/aligned_reads/")
    ids = ['data/Tox4-full/vireo/' + s + '/donor_ids.tsv' for s in samples]
  elif wildcards.group == "Tox5-full":
    samples = os.listdir("data/Tox5-full/aligned_reads/")
    ids = ['data/Tox5-full/vireo/' + s + '/donor_ids.tsv' for s in samples]
  elif wildcards.group == "Tox6-full":
    samples = os.listdir("data/Tox6-full/aligned_reads/")
    ids = ['data/Tox6-full/vireo/' + s + '/donor_ids.tsv' for s in samples]
  elif wildcards.group == "Tox7-full":
    samples = os.listdir("data/Tox7-full/aligned_reads/")
    ids = ['data/Tox7-full/vireo/' + s + '/donor_ids.tsv' for s in samples]
  return(ids)

def list_sample_summaries(wildcards):
  if wildcards.group == "Tox1":
    samples = set(glob_wildcards("data/Tox1/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples) 
    summaries = ['data/Tox1/postfiltering_read_counts/' + s + '.tsv' for s in samples]
  elif wildcards.group == "Tox1-reSEQ":
    samples = set(glob_wildcards("data/Tox1-reSEQ/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples) 
    summaries = ['data/Tox1-reSEQ/postfiltering_read_counts/' + s + '.tsv' for s in samples]
  elif wildcards.group == "Tox2":
    samples = set(glob_wildcards("data/Tox2/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples) 
    summaries = ['data/Tox2/postfiltering_read_counts/' + s + '.tsv' for s in samples]
  elif wildcards.group == "Tox3":
    samples = set(glob_wildcards("data/Tox3/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
    summaries = ['data/Tox3/postfiltering_read_counts/' + s + '.tsv' for s in samples]
  elif wildcards.group == "Tox4":
    samples = set(glob_wildcards("data/Tox4/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
    summaries = ['data/Tox4/postfiltering_read_counts/' + s + '.tsv' for s in samples]
  elif wildcards.group == "Tox5":
    samples = set(glob_wildcards("data/Tox5/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples) 
    summaries = ['data/Tox5/postfiltering_read_counts/' + s + '.tsv' for s in samples]
  elif wildcards.group == "Tox6":
    samples = set(glob_wildcards("data/Tox6/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
    summaries = ['data/Tox6/postfiltering_read_counts/' + s + '.tsv' for s in samples]
  elif wildcards.group == "Tox7":
    samples = set(glob_wildcards("data/Tox7/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
    summaries = ['data/Tox7/postfiltering_read_counts/' + s + '.tsv' for s in samples]
  return(summaries)

rule make_adata:
    resources:
        mem_mb=160000,
        time="2:00:00"
    input:
        combined_reads="data/{group}/combined_reads/combined/outs/count/filtered_feature_bc_matrix.h5",
        aggregator="data/{group}/reads_aggregator.csv",
        vireo=list_donor_ids
    output:
        outfile="data/{group}/single_cell_objects/raw.h5ad"
    params:
        prefix="data/{group}/combined_reads/combined/outs/count/filtered_feature_bc_matrix"
    conda:
        "/project2/gilad/jpopp/ebQTL/slurmy/scvi.yml"
    script: "../code/preprocessing/make_adata.py"

rule qc:
    resources:
        mem_mb=120000,
        time="1:00:00"
    input:
        full_adata="data/{group}/single_cell_objects/raw.h5ad"
    output:
        filtered_adata="data/{group}/single_cell_objects/filtered.h5ad",
        metadata="data/{group}/single_cell_objects/metadata_filtered.tsv"
    conda:
        "/project2/gilad/jpopp/ebQTL/slurmy/scvi.yml"
    script: "../code/preprocessing/qc.py"
    
rule count_reads:
    resources:
        mem_mb=120000,
        time="1:00:00"
    input:
        web_summary="data/{group}/aligned_reads/{sample}/outs/web_summary.html",
        sample_h5="data/{group}/aligned_reads/{sample}/outs/molecule_info.h5",
        metadata="data/{group}/single_cell_objects/metadata_filtered.tsv"
    output:
        postfiltering_reads="data/{group}/postfiltering_read_counts/{sample}.tsv"
    conda:
        "../envs/r-basic.yml"
    script: "../code/preprocessing/read_counts_postfilter.R"

rule count_reads_postfiltering_all:
    resources:
        mem_mb=12000,
        time="5:00"
    input:
        sample_summaries=list_sample_summaries
    output:
        full_summary="data/{group}/postfiltering_read_counts_combined/all.combined.tsv"
    shell:
        """
        cat {input} > {output}
        """

rule normalize_counts:
    resources:
        mem_mb=120000,
        time="03:00:00"
    input:
        "data/{group}/single_cell_objects/filtered.h5ad"
    output:
        "data/{group}/single_cell_objects/filtered_pflog1ppfnorm.h5ad",
        "data/{group}/single_cell_objects/filtered_pflog1ppfnorm.hvg.h5ad"
    conda: "../envs/scvi.yml"
    script:
        "../code/preprocessing/normalization.py"
        
rule scvi_embedding:
    resources:
        partition="gpu",
        mem_mb=100000,
        time="03:00:00",
        gres="gpu:2",
        nodes=2
    input:
        "data/{group}/single_cell_objects/filtered_pflog1ppfnorm.hvg.h5ad"
    output:
        "data/{group}/single_cell_objects/filtered_pflog1ppfnorm.hvg.scvi_embedding.txt",
        "data/{group}/single_cell_objects/filtered_pflog1ppfnorm.hvg.scvi_embedding.h5ad"
    params:
        "models/{group}/scvi_ldvae"
    conda: "../envs/scvi.yml"
    script:
        "../code/preprocessing/scvi_embedding.py"
        
rule umap_embedding:
    resources:
        mem_mb=100000,
        time="04:00:00"
    input:
        "data/{group}/single_cell_objects/filtered_pflog1ppfnorm.hvg.scvi_embedding.h5ad"
    output:
        "data/{group}/single_cell_objects/filtered_pflog1ppfnorm.hvg.umap_embedding.tsv",
        "data/{group}/single_cell_objects/filtered_pflog1ppfnorm.hvg.umap_embedding.h5ad",
        "data/{group}/single_cell_objects/filtered_pflog1ppfnorm.hvg.leiden.tsv"
    conda: "../envs/scvi.yml"
    script:
        "../code/preprocessing/umap_embedding.py"
        
## Below is the pipeline for the full analysis, which will hopefully replace most of the above
rule merge_adata:
    resources:
        mem_mb=750000,
        time="3:00:00",
        partition="bigmem"
    input:
        expand("data/Tox{group}-full/single_cell_objects/raw.h5ad", group=range(1, 8))
    output:
        full="data/merged_1_to_7/single_cell_objects/raw.h5ad"
    conda: "../envs/scvi.yml"
    script: "../code/preprocessing/merge_adata.py"
    
rule qc_full:
    resources:
        mem_mb=500000,
        time="30:00",
        partition="bigmem"
    input:
        full_adata="data/merged_1_to_7/single_cell_objects/raw.h5ad"
    output:
        filtered_adata="data/merged_1_to_7/single_cell_objects/filtered.h5ad",
        metadata="data/merged_1_to_7/single_cell_objects/metadata_filtered.tsv"
    conda: "../envs/scvi.yml"
    script: "../code/preprocessing/qc.py"

rule normalize_counts_full:
    resources:
        mem_mb=600000,
        time="30:00",
        partition="bigmem"
    input:
        "data/merged_1_to_7/single_cell_objects/filtered.h5ad"
    output:
        "data/merged_1_to_7/single_cell_objects/filtered_pflog1ppfnorm.h5ad",
        "data/merged_1_to_7/single_cell_objects/filtered_pflog1ppfnorm.hvg.h5ad"
    conda: "../envs/scvi.yml"
    script:
        "../code/preprocessing/normalization.py"

rule scvi_embedding_full:
    resources:
        partition="gpu",
        mem_mb=50000,
        time="05:00:00",
        gres="gpu:2",
        nodes=2
    input:
        "data/merged_1_to_7/single_cell_objects/filtered_pflog1ppfnorm.hvg.h5ad"
    output:
        "data/merged_1_to_7/single_cell_objects/filtered_pflog1ppfnorm.hvg.scvi_embedding.txt",
        "data/merged_1_to_7/single_cell_objects/filtered_pflog1ppfnorm.hvg.scvi_embedding.h5ad"
    params:
        "models/merged_1_to_7/scvi_ldvae"
    conda: "../envs/scvi.yml"
    script:
        "../code/preprocessing/scvi_embedding.py"
        
rule umap_embedding_full:
    resources:
        mem_mb=600000,
        time="03:00:00",
        partition="bigmem"
    input:
        "data/merged_1_to_7/single_cell_objects/filtered_pflog1ppfnorm.hvg.scvi_embedding.h5ad"
    output:
        "data/merged_1_to_7/single_cell_objects/filtered_pflog1ppfnorm.hvg.umap_embedding.tsv",
        "data/merged_1_to_7/single_cell_objects/filtered_pflog1ppfnorm.hvg.umap_embedding.h5ad",
        "data/merged_1_to_7/single_cell_objects/filtered_pflog1ppfnorm.hvg.leiden.tsv"
    conda: "../envs/scvi.yml"
    script:
        "../code/preprocessing/umap_embedding.py"

