# Dec 10th, 2025
# Thi Tran 
# CNV, TMB, and survival master files for Gao Feng 

# ENVIRONMENT ----
library('tidyverse')
library('ggplot2')
library('ggpubr',
        lib.loc = '/storage2/fs1/abby.green/Active/R_libraries/Bunny-Wunnies Freak Out.Bunny-Wunnies Freak Out')
working_dir <- '/storage2/fs1/abby.green/Active/Users/thi/SMC56_de_novo/2025/PCAWG'
smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
                'NSMCE4A', 'NSMCE4B', 'EID3',
                'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')

# INPUT ####
# to keep the whitelist patients only
white_list_sample <- read.csv(paste(working_dir, '/data/white_list_sample_sheet.csv', sep = ''),
                              header = TRUE)
# to keep the primary samples only 
clinical_sample <- read.delim('/storage2/fs1/abby.green/Active/RawData/TCGA_PCAWG/cBioPortal/pancan_pcawg_2020/data_clinical_sample.txt',
                              header = TRUE, sep = '\t', comment.char = '#', fill = FALSE, quote = '')
primary_sample <- filter(clinical_sample, SAMPLE_TYPE=='Primary') %>%
  inner_join(white_list_sample, by = join_by(SAMPLE_ID==icgc_specimen_id)) %>%
  filter(is.na(PLOIDY)==FALSE)
# smc56 samples
smc56_smv <- read.csv(paste(working_dir, '/data/smc56_snv.csv', sep = ''),
                      header = TRUE) %>%
  mutate(COHORT = 'SmV') %>%
  filter(Tumor_Sample_Barcode %in% primary_sample$aliquot_id)

smc56_low_cnv <- read.csv(paste(working_dir, '/data/smc56_low_cnv.csv', sep = ''),
                          header = TRUE) %>%
  mutate(COHORT = 'low CNV') %>%
  filter(Tumor_Sample_Barcode %in% primary_sample$aliquot_id)

smc56_high_cnv <- read.csv(paste(working_dir, '/data/smc56_high_cnv.csv', sep = ''),
                           header = TRUE) %>%
  mutate(COHORT = 'high CNV') %>%
  filter(Tumor_Sample_Barcode %in% primary_sample$aliquot_id)

# no introns 
smc56_smv <- filter(smc56_smv, Variant_Classification != 'Intron')
smc_ploidy <- inner_join(primary_sample, smc56_smv, 
                         by = join_by(aliquot_id == Tumor_Sample_Barcode)) %>%
  select(all_of(c('aliquot_id', 'PLOIDY', 'COHORT'))) %>%
  unique()
low_cnv_ploidy <- inner_join(primary_sample, smc56_low_cnv,
                             by = join_by(aliquot_id == Tumor_Sample_Barcode)) %>%
  select(all_of(c('aliquot_id', 'PLOIDY', 'COHORT'))) %>%
  unique()
high_cnv_ploidy <- inner_join(primary_sample, smc56_high_cnv,
                              by = join_by(aliquot_id == Tumor_Sample_Barcode)) %>%
  select(all_of(c('aliquot_id', 'PLOIDY', 'COHORT'))) %>%
  unique()
wt_ploidy <- filter(primary_sample, !aliquot_id %in% union(smc56_smv$Tumor_Sample_Barcode,
                                                           union(smc56_low_cnv$Tumor_Sample_Barcode,
                                                                 smc56_high_cnv$Tumor_Sample_Barcode))) %>%
  mutate(COHORT = 'WT') %>%
  select(all_of(c('aliquot_id', 'PLOIDY', 'COHORT'))) %>%
  unique()

df <- bind_rows(smc_ploidy, wt_ploidy,
                low_cnv_ploidy, high_cnv_ploidy) %>%
  mutate(color = with(., ifelse(COHORT %in% c('low CNV', 'high CNV'), 
                                'CNV', COHORT))) %>%
  mutate(color = forcats::fct_relevel(color, 'SmV', 'CNV', 'WT')) %>%
  select(all_of(c('PLOIDY','color','aliquot_id'))) %>%
  unique()
# output this file 
write.csv(df, paste(working_dir,'/data/no_introns_ploidy.csv',sep=''),row.names = F)
