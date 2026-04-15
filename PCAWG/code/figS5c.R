# ENVIRONMENT ####
library('tidyverse')
library('ggplot2')
library('ggpubr')
working_dir <- 'PCAWG'
theme_general <- theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_blank(),
    legend.position = 'none'
  ) 
  
# INPUT ####
# to keep the whitelist patients only
white_list_sample <- read.csv(paste(working_dir, '/data/white_list_sample_sheet.csv', sep = ''),
                              header = TRUE)
# to keep the primary samples only 
clinical_sample <- read.delim('cBioPortal/pancan_pcawg_2020/data_clinical_sample.txt',
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
# 19 March 2026 - master TMB -----
# number of tumors in this file is the number of tumors with all information available 
tmb_master <- read.csv(paste(working_dir,'/data/white_list_tmb.csv',sep=''),header=T)
primary_sample <- filter(primary_sample, aliquot_id %in% tmb_master$Tumor_Sample_Barcode)
# ANALYSIS
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
ggplot(df, aes(x=color, y=PLOIDY, color=color)) +
  geom_boxplot() +
  stat_compare_means(method = 'wilcox.test', label = 'p',
                     comparisons = list(c('CNV','WT'), c('SmV','CNV'), c('SmV','WT'))) +
  scale_color_manual(values=c('#A52A2A','#702963','#28446f'))
ggsave(paste(working_dir, '/figures/no_introns_3groups_ploidy_boxplot.eps', sep = ''),
       dpi = 320, width = 6, height = 5, units = 'in')
