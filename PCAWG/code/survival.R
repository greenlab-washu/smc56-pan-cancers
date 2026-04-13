# July 22nd, 2025
# Thi Tran 
# patients survival based on SMC5/6 alterations

# ENVIRONMENT ####
library('tidyverse')
library('ggplot2')
library('ggpubr',
        lib.loc = '/storage2/fs1/abby.green/Active/R_libraries/Bunny-Wunnies Freak Out.Bunny-Wunnies Freak Out')
library('survminer',
        lib.loc = '/storage2/fs1/abby.green/Active/R_libraries/Bunny-Wunnies Freak Out.Bunny-Wunnies Freak Out')
library('survival')

input_dir <- '/storage2/fs1/abby.green/Active/RawData/TCGA_PCAWG/cBioPortal/pancan_pcawg_2020/'
output_dir <- '/storage2/fs1/abby.green/Active/Users/thi/SMC56_de_novo/2025/PCAWG'

# INPUT ####
# load tumors information for SMC5/6 alterations
smc56_snv <- read.csv(paste(output_dir, '/data/smc56_snv.csv', sep = ''),
                      header = TRUE)
smc56_low_cnv <- read.csv(paste(output_dir, '/data/smc56_low_cnv.csv', sep = ''),
                          header = TRUE)
smc56_high_cnv <- read.csv(paste(output_dir, '/data/smc56_high_cnv.csv', sep = ''),
                           header = TRUE)
raw_clinical_sample <- read.delim(paste(input_dir, 'data_clinical_sample.txt', sep = ''),
                                  header = TRUE, sep = '\t', comment.char = '#', fill = FALSE, quote = '')
raw_clinical_patients <- read.delim(paste(input_dir, 'data_clinical_patient.txt', sep = ''),
                                  header = TRUE, sep = '\t', comment.char = '#', fill = FALSE, quote = '')
white_list_sample <- read.csv(paste(output_dir, '/data/white_list_sample_sheet.csv', sep = ''),
                         header = TRUE)

# CLEANUP and ANALYSIS ####
snv_sample <- filter(sample_sheet, aliquot_id %in% smc56_snv$Tumor_Sample_Barcode) %>%
  unique()
cnv_sample <- filter(sample_sheet, aliquot_id %in% union(
  smc56_low_cnv$Tumor_Sample_Barcode, smc56_high_cnv$Tumor_Sample_Barcode
)) %>%
  unique()
wt_sample <- sample_sheet %>%
  filter(!(icgc_donor_id %in% union(
  snv_sample$icgc_donor_id, cnv_sample$icgc_donor_id
)))

snv_patients <- select(raw_clinical_patients, all_of(c('PATIENT_ID', 'OS_STATUS', 'OS_MONTHS'))) %>%
  unique() %>%
  filter(PATIENT_ID %in% snv_sample$icgc_donor_id) %>%
  mutate(COHORT = 'SmV')
cnv_patients <- select(raw_clinical_patients, all_of(c('PATIENT_ID', 'OS_STATUS', 'OS_MONTHS'))) %>%
  unique() %>%
  filter(PATIENT_ID %in% cnv_sample$icgc_donor_id) %>%
  mutate(COHORT = 'CNV')
wt_patients <- select(raw_clinical_patients, all_of(c('PATIENT_ID', 'OS_STATUS', 'OS_MONTHS'))) %>%
  unique() %>%
  filter(PATIENT_ID %in% wt_sample$icgc_donor_id) %>%
  mutate(COHORT = 'WT')

df <- rbind(snv_patients, cnv_patients, wt_patients) %>% 
  filter(is.na(OS_MONTHS)==FALSE) %>%
  separate(OS_STATUS, into = c('OS_STATUS_NUM', 'OS_STATUS_STATE'), sep = ':') %>%
  mutate(OS_STATUS_NUM = as.numeric(OS_STATUS_NUM))
# construct Kaplan-Meier survival table for each cohort
fit <- df %>%
  surv_fit(Surv(OS_MONTHS, OS_STATUS_NUM) ~ COHORT, data = .)

# GRAPH ####
# control the aesthetic for the graph
theme_general <- theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

ggsurvplot(fit, data = df, pval = TRUE,
           ggtheme = theme_general,
           legend = 'right', legend.title = '',
           xlab = 'Time (month)')
ggsave(paste(output_dir, '/figures/survival_3groups_v1.eps', sep = ''),
       dpi = 320, width = 5, height = 4.5, units = 'in')

# PRIMARY only - no introns ####
clinical_sample <- read.delim('/storage2/fs1/abby.green/Active/RawData/TCGA_PCAWG/cBioPortal/pancan_pcawg_2020/data_clinical_sample.txt',
                              header = TRUE, sep = '\t', comment.char = '#', fill = FALSE, quote = '')
primary_sample <- filter(clinical_sample, SAMPLE_TYPE=='Primary') %>%
  inner_join(white_list_sample, by = join_by(SAMPLE_ID==icgc_specimen_id)) %>%
  filter(is.na(PLOIDY)==FALSE)

smc56_snv <- filter(smc56_snv, Variant_Classification != 'Intron')
primary_snv_sample <- filter(primary_sample, aliquot_id %in% smc56_snv$Tumor_Sample_Barcode) %>%
  unique()
primary_cnv_sample <- filter(primary_sample, aliquot_id %in% union(
  smc56_low_cnv$Tumor_Sample_Barcode, smc56_high_cnv$Tumor_Sample_Barcode
)) %>%
  unique()
primary_wt_sample <- filter(primary_sample, !(icgc_donor_id %in% union(
    primary_snv_sample$icgc_donor_id, primary_cnv_sample$icgc_donor_id
  ))) %>%
  unique()

snv_patients <- select(raw_clinical_patients, all_of(c('PATIENT_ID', 'OS_STATUS', 'OS_MONTHS'))) %>%
  unique() %>%
  filter(PATIENT_ID %in% primary_snv_sample$icgc_donor_id) %>%
  mutate(COHORT = 'SmV')
cnv_patients <- select(raw_clinical_patients, all_of(c('PATIENT_ID', 'OS_STATUS', 'OS_MONTHS'))) %>%
  unique() %>%
  filter(PATIENT_ID %in% primary_cnv_sample$icgc_donor_id) %>%
  mutate(COHORT = 'CNV')
wt_patients <- select(raw_clinical_patients, all_of(c('PATIENT_ID', 'OS_STATUS', 'OS_MONTHS'))) %>%
  unique() %>%
  filter(PATIENT_ID %in% primary_wt_sample$icgc_donor_id) %>%
  mutate(COHORT = 'WT')

df <- rbind(snv_patients, cnv_patients, wt_patients) %>% 
  filter(is.na(OS_MONTHS)==FALSE) %>%
  separate(OS_STATUS, into = c('OS_STATUS_NUM', 'OS_STATUS_STATE'), sep = ':') %>%
  mutate(OS_STATUS_NUM = as.numeric(OS_STATUS_NUM))

fit <- df %>%
  surv_fit(Surv(OS_MONTHS, OS_STATUS_NUM) ~ COHORT, data = .)
pairwise_survdiff(Surv(OS_MONTHS, OS_STATUS_NUM) ~ COHORT, data = df,
                  p.adjust.method = 'none')

ggsurvplot(fit, data = df, pval = TRUE,
           ggtheme = theme_general,
           legend = 'right', legend.title = '',
           xlab = 'Time (month)')
ggsave(paste(output_dir, '/figures/no_introns_3groups_survival.eps', sep = ''),
       dpi= 320, width = 6, height = 5)
