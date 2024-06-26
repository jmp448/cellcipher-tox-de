rule germlayer_tox_de_celltype:
    resources:
        mem_mb=60000,
        time="1:00:00"
    input:
        "data/training/tox_de_training_data.RData"
    output:
        "results/germlayer/tox_de.{celltype}.rds"
    conda:
        "../slurmy/r-de.yml"
    script: "../code/tox_de/tox_de.germlayer.R"

rule germlayer_sitaxsentan_de_celltype:
    resources:
        mem_mb=60000,
        time="1:00:00"
    input:
        "data/training/tox_de_training_data.RData"
    output:
        "results/germlayer/sitaxsentan_de.{celltype}.rds"
    conda:
        "../slurmy/r-de.yml"
    script: "../code/tox_de/sitaxsentan_de.germlayer.R"

rule germlayer_tolcapone_de_celltype:
    resources:
        mem_mb=60000,
        time="1:00:00"
    input:
        "data/training/tox_de_training_data.RData"
    output:
        "results/germlayer/tolcapone_de.{celltype}.rds"
    conda:
        "../slurmy/r-de.yml"
    script: "../code/tox_de/tolcapone_de.germlayer.R"

rule germlayer_nefazodone_de_celltype:
    resources:
        mem_mb=60000,
        time="1:00:00"
    input:
        "data/training/tox_de_training_data.RData"
    output:
        "results/germlayer/nefazodone_de.{celltype}.rds"
    conda:
        "../slurmy/r-de.yml"
    script: "../code/tox_de/nefazodone_de.germlayer.R"

rule germlayer_troglitazone_de_celltype:
    resources:
        mem_mb=60000,
        time="1:00:00"
    input:
        "data/training/tox_de_training_data.RData"
    output:
        "results/germlayer/troglitazone_de.{celltype}.rds"
    conda:
        "../slurmy/r-de.yml"
    script: "../code/tox_de/troglitazone_de.germlayer.R"

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

