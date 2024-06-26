library(DropletUtils)
library(tidyverse)
library(data.table)

# molecule_h5 <- '/project2/gilad/jpopp/cellcipher/data_071023/aligned_reads/YG-KR-10X-13s-Tox-01/outs/molecule_info.h5'
# sample <- 'YG-KR-10X-13s-Tox-01'

molecule_h5 <- snakemake@input[['sample_h5']]
metadata_loc <- snakemake@input[['metadata']]
web_summary_file <- snakemake@input[['web_summary']]
sample <- snakemake@wildcards[['sample']]
postfiltering_reads_loc <- snakemake@output[['postfiltering_reads']]

# Load web summary data
web_summary <- readLines(web_summary_file)
total_reads_line <- grep('"Number of Reads"', web_summary, value = TRUE)
total_reads_match <- str_match(total_reads_line, '"Number of Reads", "(\\d+(,\\d+)*)"')
total_reads <- as.numeric(gsub(",", "", total_reads_match[1, 2]))

# Load molecule info
molecule_info <- read10xMolInfo(molecule_h5)
sample_metadata <- read_tsv(metadata_loc) %>%
  filter(sample_id==sample) %>%
  mutate(cell=str_extract(cell, "[^-]+"))

# Count the total number of aligned reads in the sample
h5_table <- as.data.table(molecule_info$data@listData)
total_aligned_reads <- sum(h5_table$reads)

# Count the reads that make it past filtering
h5_filter <- h5_table[h5_table$cell %in% sample_metadata$cell]
total_postfiltering_reads <- sum(h5_filter$reads)

summary <- tibble("sample"=sample, 
                  "total_reads"=total_reads,
                  "total_aligned_reads"=total_aligned_reads,
                  "total_postfiltering_reads"=total_postfiltering_reads) %>%
  write_tsv(postfiltering_reads_loc)
