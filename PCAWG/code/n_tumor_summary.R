# ENVIRONMENT ####
library('tidyverse')
library('ggplot2')
input_dir <- '/storage2/fs1/abby.green/Active/RawData/TCGA_PCAWG/PCAWG/'
output_dir <- '/storage2/fs1/abby.green/Active/Users/thi/SMC56_de_novo/2025/PCAWG'

# INPUT ####
raw_snv <- read.delim(paste(input_dir, '/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf',
                            sep = ''), header = TRUE, sep = '\t', comment.char = '#',
                      fill = FALSE, quote = '')
raw_cnv <- read.delim(paste(input_dir, '/consensus_cnv/gene_level_calls/all_samples.consensus_CN.by_gene.170214.txt',
                            sep = ''), header = TRUE, sep = '\t', fill = FALSE, comment.char = '#', quote = '')
colnames(raw_cnv) <- stringr::str_replace_all(colnames(raw_cnv), '\\.', '_') %>%
  stringr::str_replace('^X', '')
icgc_clinical_1 <- read.delim(paste(input_dir, '/pcawg_sample_sheet.tsv', sep = ''),
                              header = TRUE, sep = '\t', comment.char = '#',
                              fill = FALSE, quote = '') %>%
  separate_wider_delim(donor_unique_id, delim = '::', names = c('project', 'donor_unique_id'))
icgc_clinical_2 <- read.delim(paste(input_dir, '/specimen.all_projects.tsv', sep = ''),
                              header = TRUE, sep = '\t', comment.char = '#',
                              fill = FALSE, quote = '')
icgc_clinical_3 <- read.csv(paste(input_dir, 'pcawg_donor_clinical_August2016_v9.csv', sep = ''), 
                            header = TRUE)

cbioportal_clinical_sample <- read.delim('/storage2/fs1/abby.green/Active/RawData/TCGA_PCAWG/cBioPortal/pancan_pcawg_2020/data_clinical_sample.txt',
                                         header = TRUE, sep = '\t', comment.char = '#', fill = FALSE, quote = '')
cbioportal_clinical_patient <- read.delim('/storage2/fs1/abby.green/Active/RawData/TCGA_PCAWG/cBioPortal/pancan_pcawg_2020/data_clinical_patient.txt',
                                          header = TRUE, sep = '\t', comment.char = '#', fill = FALSE, quote = '')
# keep the whitelist patients only 
white_sample_sheet <- read.csv(paste(output_dir, '/data/white_list_sample_sheet.csv',sep=''),
                               header=T)
smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
                'NSMCE4A', 'NSMCE4B', 'EID3',
                'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')

patient_data <- filter(icgc_clinical_2, PATIENT_ID %in% white_sample_sheet$icgc_donor_id)
# Tumors with SMC56 SNVs and indels 
smc56_snv <- filter(raw_snv, 
                    Hugo_Symbol %in% smc56_gene & Variant_Classification != 'Intron')
# Tumors with SMC56 CNVs 
smc56_low_cnv <- filter(raw_cnv, Gene_Symbol %in% smc56_gene) %>% # filter based on SMC5/6 copy number
  select(-any_of(c('Locus_ID', 'Cytoband'))) %>% # adjust the dataframe so I can filter based on gene level copy number
  pivot_longer(cols = -Gene_Symbol, names_to = 'Tumor_Sample_Barcode', values_to = 'copy_number') %>%
  filter(copy_number < 2 | copy_number == 'NaN') %>%
  mutate_at(c('Tumor_Sample_Barcode'), stringr::str_remove, '^X') %>% # adjust the Tumor_Sample_Barcode column for matching later
  mutate_at(c('Tumor_Sample_Barcode'), stringr::str_replace_all, '_', '-')
# Tumors with high SMC5/6 CNVs 
smc56_high_cnv <- filter(raw_cnv, Gene_Symbol %in% smc56_gene) %>% # filter based on SMC5/6 copy number
  select(-any_of(c('Locus_ID', 'Cytoband'))) %>% # adjust the dataframe so I can filter based on gene level copy number
  pivot_longer(cols = -Gene_Symbol, names_to = 'Tumor_Sample_Barcode', values_to = 'copy_number') %>%
  filter(copy_number != 'NaN' & copy_number > 2) %>%
  mutate_at(c('Tumor_Sample_Barcode'), stringr::str_remove, '^X') %>% # adjust the Tumor_Sample_Barcode column for matching later
  mutate_at(c('Tumor_Sample_Barcode'), stringr::str_replace_all, '_', '-')
smc56_cnv <- filter(white_sample_sheet, 
                    aliquot_id %in% union(smc56_high_cnv$Tumor_Sample_Barcode, smc56_low_cnv$Tumor_Sample_Barcode))
