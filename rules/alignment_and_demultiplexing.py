#Snakefile 
# TODO update so that tox2_samples / tox3_samples isn't like manually updated
# TODO update so that the donors included aren't just the same 5

import glob
import os
import sys

# Functions
# def list_samples(wildcards):
#   if wildcards.group == "Tox1-full":
#     samples = set(glob_wildcards("data/Tox1-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples) 
#   elif wildcards.group == "Tox2-full":
#     samples = set(glob_wildcards("data/Tox2-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples) 
#   elif wildcards.group == "Tox3-full":
#     samples = set(glob_wildcards("data/Tox3-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
#   elif wildcards.group == "Tox4-full":
#     samples = set(glob_wildcards("data/Tox4-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
#   elif wildcards.group == "Tox5-full":
#     samples = set(glob_wildcards("data/Tox5-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples) 
#   elif wildcards.group == "Tox6-full":
#     samples = set(glob_wildcards("data/Tox6-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
#   elif wildcards.group == "Tox7-full":
#     samples = set(glob_wildcards("data/Tox7-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
#   return(samples)

def list_samples(wildcards):
  if wildcards.group == "Tox1-full":
    samples = os.listdir("data/Tox1-full/aligned_reads/") # note we don't want to include 'unassigned', which is in fastqs
  elif wildcards.group == "Tox2-full":
    samples = set(glob_wildcards("data/Tox2-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
  elif wildcards.group == "Tox3-full":
    samples = set(glob_wildcards("data/Tox3-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
  elif wildcards.group == "Tox4-full":
    samples = set(glob_wildcards("data/Tox4-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
  elif wildcards.group == "Tox5-full":
    samples = set(glob_wildcards("data/Tox5-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
  elif wildcards.group == "Tox6-full":
    samples = set(glob_wildcards("data/Tox6-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
  elif wildcards.group == "Tox7-full":
    samples = set(glob_wildcards("data/Tox7-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
  return(filter_to_unique(samples))

def list_files_to_transfer(wildcards):
  if wildcards.group == "Tox1":
    fastqs_orig = list(glob_wildcards("data/Tox1/fastqs/{fastas}.fastq.gz").fastas)
    fastqs_orig = [f for f in fastqs_orig if not f.startswith("Undetermined")]
    fastqs_new = ["data/Tox1-full/fastqs/" + f + ".fastq.gz" for f in fastqs_orig]
  elif wildcards.group == "Tox2":
    fastqs_orig = list(glob_wildcards("data/Tox2/fastqs/{fastas}.fastq.gz").fastas)
    fastqs_new = ["data/Tox2-full/fastqs/" + f + ".fastq.gz" for f in fastqs_orig]
  elif wildcards.group == "Tox3":
    fastqs_orig = list(glob_wildcards("data/Tox3/fastqs/{fastas}.fastq.gz").fastas)
    fastqs_new = ["data/Tox3-full/fastqs/" + f + ".fastq.gz" for f in fastqs_orig]
  elif wildcards.group == "Tox4":
    fastqs_orig = list(glob_wildcards("data/Tox4/fastqs/{fastas}.fastq.gz").fastas)
    fastqs_new = ["data/Tox4-full/fastqs/" + f + ".fastq.gz" for f in fastqs_orig]
  elif wildcards.group == "Tox5":
    fastqs_orig = list(glob_wildcards("data/Tox5/fastqs/{fastas}.fastq.gz").fastas)
    fastqs_new = ["data/Tox5-full/fastqs/" + f + ".fastq.gz" for f in fastqs_orig]
  elif wildcards.group == "Tox6":
    fastqs_orig = list(glob_wildcards("data/Tox6/fastqs/{fastas}.fastq.gz").fastas)
    fastqs_new = ["data/Tox6-full/fastqs/" + f + ".fastq.gz" for f in fastqs_orig]
  elif wildcards.group == "Tox7":
    fastqs_orig = list(glob_wildcards("data/Tox7/fastqs/{fastas}.fastq.gz").fastas)
    fastqs_new = ["data/Tox7-full/fastqs/" + f + ".fastq.gz" for f in fastqs_orig]
  return(fastqs_new)

# def list_molecule_h5_files(wildcards):
#   if wildcards.group == "Tox1-full":
#     samples = set(glob_wildcards("data/Tox1-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples) 
#     h5s = ['data/Tox1-full/aligned_reads/' + s + '/outs/molecule_info.h5' for s in samples]
#   elif wildcards.group == "Tox2-full":
#     samples = set(glob_wildcards("data/Tox2-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples) 
#     h5s = ['data/Tox2-full/aligned_reads/' + s + '/outs/molecule_info.h5' for s in samples]
#   elif wildcards.group == "Tox3-full":
#     samples = set(glob_wildcards("data/Tox3-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
#     h5s = ['data/Tox3-full/aligned_reads/' + s + '/outs/molecule_info.h5' for s in samples]
#   elif wildcards.group == "Tox4-full":
#     samples = set(glob_wildcards("data/Tox4-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
#     h5s = ['data/Tox4-full/aligned_reads/' + s + '/outs/molecule_info.h5' for s in samples]
#   elif wildcards.group == "Tox5-full":
#     samples = set(glob_wildcards("data/Tox5-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
#     h5s = ['data/Tox5-full/aligned_reads/' + s + '/outs/molecule_info.h5' for s in samples]
#   elif wildcards.group == "Tox6-full":
#     samples = set(glob_wildcards("data/Tox6-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
#     h5s = ['data/Tox6-full/aligned_reads/' + s + '/outs/molecule_info.h5' for s in samples]
#   elif wildcards.group == "Tox7-full":
#     samples = set(glob_wildcards("data/Tox7-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
#     h5s = ['data/Tox7-full/aligned_reads/' + s + '/outs/molecule_info.h5' for s in samples]
#   return(h5s)

def list_samples_fcmerge(wildcards):
  if wildcards.group == "Tox1-full":
    samples = set(glob_wildcards("data/Tox1-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
  elif wildcards.group == "Tox2-full":
    samples = set(glob_wildcards("data/Tox2-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
  elif wildcards.group == "Tox3-full":
    samples = set(glob_wildcards("data/Tox3-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
  elif wildcards.group == "Tox4-full":
    samples = set(glob_wildcards("data/Tox4-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
  elif wildcards.group == "Tox5-full":
    samples = set(glob_wildcards("data/Tox5-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
  elif wildcards.group == "Tox6-full":
    samples = set(glob_wildcards("data/Tox6-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
  elif wildcards.group == "Tox7-full":
    samples = set(glob_wildcards("data/Tox7-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
  sample_fcs = []
  for s in samples:
      sample_full = s.split('-')
      sample_key = '-'.join(sample_full[:7])
      if sample_key == wildcards.sample:
        sample_fcs.append(s)
  return([','.join(sample_fcs)])

def filter_to_unique(samples):
  unique_samples = []
  for s in samples:
      sample_full = s.split('-')
      sample_key = '-'.join(sample_full[:7])
      if sample_key not in unique_samples:
        unique_samples.append(sample_key)
  return(unique_samples)
  
def list_molecule_h5_files(wildcards):
  if wildcards.group == "Tox1-full":
    samples = os.listdir("data/Tox1-full/aligned_reads/")
    h5s = ['data/Tox1-full/aligned_reads/' + s + '/outs/molecule_info.h5' for s in filter_to_unique(samples)]
  elif wildcards.group == "Tox2-full":
    samples = set(glob_wildcards("data/Tox2-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
    h5s = ['data/Tox2-full/aligned_reads/' + s + '/outs/molecule_info.h5' for s in filter_to_unique(samples)]
  elif wildcards.group == "Tox3-full":
    samples = set(glob_wildcards("data/Tox3-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
    h5s = ['data/Tox3-full/aligned_reads/' + s + '/outs/molecule_info.h5' for s in filter_to_unique(samples)]
  elif wildcards.group == "Tox4-full":
    samples = set(glob_wildcards("data/Tox4-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
    h5s = ['data/Tox4-full/aligned_reads/' + s + '/outs/molecule_info.h5' for s in filter_to_unique(samples)]
  elif wildcards.group == "Tox5-full":
    samples = set(glob_wildcards("data/Tox5-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
    h5s = ['data/Tox5-full/aligned_reads/' + s + '/outs/molecule_info.h5' for s in filter_to_unique(samples)]
  elif wildcards.group == "Tox6-full":
    samples = set(glob_wildcards("data/Tox6-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
    h5s = ['data/Tox6-full/aligned_reads/' + s + '/outs/molecule_info.h5' for s in filter_to_unique(samples)]
  elif wildcards.group == "Tox7-full":
    samples = set(glob_wildcards("data/Tox7-full/fastqs/{samples}_S{order}_L{lane}_R1_{chunk}.fastq.gz").samples)
    h5s = ['data/Tox7-full/aligned_reads/' + s + '/outs/molecule_info.h5' for s in filter_to_unique(samples)]
  return(h5s)

# def list_donors(wildcards):
#   if wildcards.group == "Tox1":
#     donors=['NA18519,NA19210,NA18912,NA19204,NA19093']
#   elif wildcards.group == "Tox2":
#     donors=['NA18519,NA19210,NA18912,NA19204,NA19093']
#   elif wildcards.group == "Tox3":
#     donors=['NA18519,NA19210,NA18912,NA19204,NA19093']
#   elif wildcards.gropu == "Tox4":
#     donors=['NA18519,NA19210,NA18912,NA19204,NA19093']

# Configuration ----------------------------------------------------------------
cellranger="/project2/gilad/kenneth/software/cellranger-5.0.1/cellranger" #cellranger executable
cellranger_aggr="/project2/gilad/kenneth/software/cellranger-6.1.2/cellranger"
#############################################################################
# Alignment
#############################################################################

rule merge_fastqs_reseq:
    input:
        orig_fastq="data/{group}/fastqs/{f}.fastq.gz"
    output:
        copy_fastq="data/{group}-full/fastqs/{f}.fastq.gz"
    run:
        if os.path.isfile(str(output.copy_fastq)):
            sys.exit("may need to overwrite {input}")
        else:
            shell("cp {input} {output}")

rule merge_fastqs_reseq_all:
    input:
        fastqs=list_files_to_transfer
    output:
        "temp/{group}-full.merged"
    shell:
      """echo done > {output}"""
            
rule alignment:
    resources:
        partition="caslake",
        mem_mb=50000,
        time="10:00:00",
        ntasks_per_node=28
    output:
        matrix = "data/{group}/aligned_reads/{sample}/outs/filtered_feature_bc_matrix.h5",
        barcodes = "data/{group}/aligned_reads/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        bam = "data/{group}/aligned_reads/{sample}/outs/possorted_genome_bam.bam",
        molecule_h5="data/{group}/aligned_reads/{sample}/outs/molecule_info.h5"
    params:
        fastq = "../fastqs",
        genomeidx = "/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A",
        fc_merged_samples=list_samples_fcmerge
    benchmark:
        "benchmarks/cellranger.{group}.{sample}.txt"
    shell:
        """
        [ ! -f data/{wildcards.group}/aligned_reads/{wildcards.sample} ] && rm -r data/{wildcards.group}/aligned_reads/{wildcards.sample}
        mkdir -p data/{wildcards.group}/aligned_reads/{wildcards.sample} && cd data/{wildcards.group}/aligned_reads/
        {cellranger} count --id={wildcards.sample} \
                 --transcriptome={params.genomeidx} \
                 --fastqs={params.fastq} \
                 --sample={params.fc_merged_samples} \
                 --chemistry=SC3Pv3 \
                 --expect-cells=10000 \
                 --include-introns=true \
                 --localcores=28 \
                 --nosecondary
        """

rule prep_csv_to_aggregate:
    input:
        fastqs=list_molecule_h5_files
        #fastqs=expand("data/{{group}}/aligned_reads/{sample}/outs/molecule_info.h5", sample=os.listdir("data/{group}/aligned_reads/"))
    output:
        "data/{group}/reads_aggregator.csv"
    params:
        samples=list_samples
    run:
        with open(output[0], "w") as f:
            f.write(f"sample_id,molecule_h5\n")
            for s, fname in zip(params.samples, input.fastqs):
                fname_redirect = fname.replace(f"data/{wildcards.group}/", "../")
                f.write(f"{s},{fname_redirect}\n")
        
rule aggregate:
    resources:
        mem_mb=100000,
        time="06:00:00"
    input:
        reads=list_molecule_h5_files,
        aggregator="data/{group}/reads_aggregator.csv"
    output:
        "data/{group}/combined_reads/combined/outs/count/filtered_feature_bc_matrix.h5"
    shell:
        """
        [ ! -f data/{wildcards.group}/combined_reads/combined ] && rm -r data/{wildcards.group}/combined_reads/combined
        mkdir -p data/{wildcards.group}/combined_reads/combined && cd data/{wildcards.group}/combined_reads/
        {cellranger_aggr} aggr --id=combined \
                 --csv=../../../{input.aggregator} \
                 --normalize=none \
                 --nosecondary \
                 --disable-ui
        """

#############################################################################
# Demultiplexing
#############################################################################

# subset vcfs
rule subset_vcf_for_vireo:
    resources:
        mem_mb=15000
    input:
        vcf="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz"
    output:
        unsorted="data/genotypes/tox_donors.vireo.unsorted.vcf",
        sorted="data/genotypes/tox_donors.vireo.vcf",
        bed="data/genotypes/tox_donors.vireo.vcf.bed"
    params:
        donors="NA18519,NA19210,NA18912,NA19204,NA19093",
        compressed ="data/genotypes/tox_donors.vireo.vcf.gz"
    shell:
        """
        module load bcftools bedtools
        source code/alignment/filter_vcf_file_for_popscle.sh
        subset_samples_from_vcf {params.donors} {input.vcf} \
          | filter_out_mutations_missing_genotype_for_one_or_more_samples \
          | filter_out_mutations_homozygous_reference_in_all_samples \
          | filter_out_mutations_homozygous_in_all_samples \
          > {output.unsorted}
        bcftools sort -Oz {output.unsorted} -o {params.compressed}
        gunzip {params.compressed}
        bedtools merge -i {output.sorted} > {output.bed}
        """
        
rule run_cellsnp:
    resources:
        mem_mb=10000,
        time="10:00:00",
        ntasks_per_node=8
    input:
         vcf = "data/genotypes/tox_donors.vireo.vcf",
         bam = "data/{group}/aligned_reads/{sample}/outs/possorted_genome_bam.bam",
         barcodes = "data/{group}/aligned_reads/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
    output:
        "data/{group}/cellsnp/{sample}/cellSNP.base.vcf.gz"
    params:
         barcodes = "temp/{sample}_barcodes.txt",
         outdir = "data/{group}/cellsnp/{sample}/"
    conda:
        "../envs/vireo.yaml"
    shell:
        """
        mkdir -p data/{wildcards.group}/cellsnp/{wildcards.sample}
        gunzip < {input.barcodes} > {params.barcodes}
        cellsnp-lite -s {input.bam} -b {params.barcodes} -O {params.outdir} -R {input.vcf} -p 8 --minMAF 0.1 --minCOUNT 20 --gzip
        """
        
rule run_vireo:
    resources:
        mem_mb=50000,
        time="36:00:00",
        ntasks_per_node=4
    input:
        vcf = "data/genotypes/tox_donors.vireo.vcf",
        cellsnp = "data/{group}/cellsnp/{sample}/cellSNP.base.vcf.gz"
    output:
         "data/{group}/vireo/{sample}/donor_ids.tsv"
    params:
         vireo = "data/{group}/vireo/{sample}",
         cellsnp = "data/{group}/cellsnp/{sample}",
         nind = 5
    conda:
        "../envs/vireo.yaml"
    shell:
        """
        mkdir -p data/{wildcards.group}/vireo/{wildcards.sample}
        vireo -c {params.cellsnp} -d {input.vcf} -o {params.vireo} -t GT --nproc=4 -N {params.nind} --noPlot
        """
