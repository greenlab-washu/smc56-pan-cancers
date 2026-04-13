# June 26th, 2025
# Thi Tran 
# Generate 96-trinucleotide context and SBSsig for PCAWG
# from the individual vcf.gz files 

# ENVIRONMENT ####
#library('tidyverse')
#library('ggplot2')
library('MutationalPatterns')
library('BSgenome')
library('dplyr')
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#BiocManager::install(ref_genome)
library(ref_genome, character.only = TRUE)

input_dir <- '/storage2/fs1/abby.green/Active/RawData/TCGA_PCAWG/PCAWG/consensus_snv_indel/snv_indel_icgc/snv_mnv'
output_dir <- '/storage2/fs1/abby.green/Active/Users/thi/SMC56_de_novo/2025/PCAWG'

# INPUT ####
vcf_files <- list.files(
  path = input_dir, 
  full.names = TRUE, pattern = 'gz$'
)
sample_names <- stringr::str_extract(vcf_files, '(?<=snv_mnv/).+(?=\\.consensus)')
grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)

# 96 trinucleotide context and output it ####
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome) %>%
  as.data.frame()
write.csv(mut_mat, paste(output_dir, '/data/indv_96context.csv', sep = ''),
          row.names = FALSE)
