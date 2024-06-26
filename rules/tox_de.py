rule tox_de_glmgampoi:
    resources:
        mem_mb=60000,
        time="8:00:00",
        nodes=1,
        ntasks_per_node=1
    input:
        "data/training/tox_de_training_data.RData"
    output:
        "results/tox_de.glmgampoi.rds"
    conda:
        "../slurmy/r-de.yml"
    script: "../code/tox_de/tox_DE_glmgampoi.R"
    
rule tox_de_lrt:
    resources:
        mem_mb=60000,
        time="1:00:00"
    input:
        "data/training/tox_de_training_data.RData"
    output:
        "results/tox_de.lrt.rds"
    conda:
        "../slurmy/r-de.yml"
    script: "../code/tox_de/tox_DE_lrt.R"

rule tox_de_lrt_celltype:
    resources:
        mem_mb=60000,
        time="1:00:00"
    input:
        "data/training/tox_de_training_data.RData"
    output:
        "results/tox_de.lrt.{celltype}.rds"
    conda:
        "../slurmy/r-de.yml"
    script: "../code/tox_de/tox_DE_lrt_celltype.R"

rule tox_de_mash:
    resources:
        mem_mb=60000,
        time="1:00:00"
    input:
        "results/tox_de.mash_inputs.RData"
    output:
        "results/tox_de.mash_outputs.RData"
    conda:
        "../slurmy/r-mashr.yml" ##NOTE: mash installation failed when first creating this env, but the conda installation said it succeeded, so I opened the conda env, started R, and installed mash from CRAN
    script: "../code/tox_de/tox_DE_mash.R"
    
rule tox_de_lrt_donor_celltype_intx:
    resources:
        mem_mb=60000,
        time="1:00:00"
    input:
        "data/training/tox_de_training_data.RData"
    output:
        "results/tox_de.lrt.donor_celltype_intx.rds"
    conda:
        "../slurmy/r-de.yml"
    script: "../code/tox_de/tox_DE_lrt_donor_celltype_intx.R"

rule tox_de_lrt_donor_celltype_intx_permutations:
    resources:
        mem_mb=60000,
        time="30:00"
    input:
        "data/training/tox_de_training_data.RData"
    output:
        "results/tox_de.lrt.donor_celltype_intx.{seed}.tsv"
    conda:
        "../slurmy/r-de.yml"
    script: "../code/tox_de/tox_DE_lrt_donor_celltype_intx.R"
