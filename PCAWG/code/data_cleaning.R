# define tumors with SMC5/6 genetic alterations in PCAWG dataset

# ENVIRONMENT ####
library('tidyverse')
library('ggplot2')
input_dir <- 'PCAWG'
output_dir <- 'PCAWG'

# INPUT ####
raw_snv <- read.delim(paste(input_dir, '/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf',
                            sep = ''), header = TRUE, sep = '\t', comment.char = '#',
                      fill = FALSE, quote = '')
raw_cnv <- read.delim(paste(input_dir, '/consensus_cnv/gene_level_calls/all_samples.consensus_CN.by_gene.170214.txt',
                            sep = ''), header = TRUE, sep = '\t', fill = FALSE, comment.char = '#', quote = '')
colnames(raw_cnv) <- stringr::str_replace_all(colnames(raw_cnv), '\\.', '-') %>%
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

cbioportal_clinical_sample <- read.delim('cBioPortal/pancan_pcawg_2020/data_clinical_sample.txt',
                           header = TRUE, sep = '\t', comment.char = '#', fill = FALSE, quote = '')
cbioportal_clinical_patient <- read.delim('cBioPortal/pancan_pcawg_2020/data_clinical_patient.txt',
                           header = TRUE, sep = '\t', comment.char = '#', fill = FALSE, quote = '')

smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
          'NSMCE4A', 'NSMCE4B', 'EID3',
          'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')
# ANALYSIS ####
## define tumors that belong in different cohorts ----
# keep the whitelist patients only 
white_sample_sheet <- filter(icgc_clinical_1, donor_wgs_exclusion_white_gray == 'Whitelist')
white_sample_sheet %>% select(all_of(c('aliquot_id', 'icgc_donor_id',
                                       'icgc_specimen_id', 'icgc_sample_id'))) %>%
  write.csv(paste(output_dir, '/data/white_list_sample_sheet.csv', sep = ''),
            row.names = FALSE)
# whitelist patients with clinical data available
patient_data <- filter(icgc_clinical_3, 
                       icgc_donor_id %in% white_sample_sheet$icgc_donor_id &
                         is.na(donor_survival_time)==F) 
# make sure that they also have CNV and SNV information as well
df <- filter(white_sample_sheet,
             aliquot_id %in% intersect(colnames(raw_cnv),unique(raw_snv$Tumor_Sample_Barcode)) &
               icgc_donor_id %in% patient_data$icgc_donor_id)
length(unique(df$icgc_donor_id)) 
# Tumors with SMC56 SNVs and indels 
smc56_snv <- filter(raw_snv, Hugo_Symbol %in% smc56_gene)
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
## calculate TMB for all tumors ----
master_tmb <- select(raw_snv, all_of(c('Tumor_Sample_Barcode', 'Hugo_Symbol',
                                       'Chromosome', 'Start_position', 'End_position',
                                       'Genome_Change'))) %>%
  unique() %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(abs_mut_number=n()) %>%
  mutate(TMB = abs_mut_number/2800) 
## recurrent SMC5/6 mutations ----
smc56_recurrent <- smc56_snv %>%
  group_by(Hugo_Symbol, Tumor_Sample_Barcode, ) %>%
# OUTPUT ####
## define tumors IDs in each cohort 
write.csv(smc56_snv, paste(output_dir, '/data/smc56_snv.csv', sep = ''),
          row.names = FALSE)
write.csv(smc56_low_cnv, paste(output_dir, '/data/smc56_low_cnv.csv', sep = ''),
          row.names = FALSE)
write.csv(smc56_high_cnv, paste(output_dir, '/data/smc56_high_cnv.csv',sep = ''),
          row.names = FALSE)

## TMB of all tumors 
write.csv(master_tmb, paste(output_dir, '/data/white_list_tmb.csv', sep = ''),
          row.names = FALSE)