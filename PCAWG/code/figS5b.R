# TMB of 3 groups: SmV, CNV, WT

# ENVIRONMENT ####
library('tidyverse')
library('ggplot2')
library('ggpubr')
working_dir <- 'PCAWG'
smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
                'NSMCE4A', 'NSMCE4B', 'EID3',
                'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')

# INPUT ####
# to keep the whitelist patients only
white_list_sample <- read.csv(paste(working_dir, '/data/white_list_sample_sheet.csv', sep = ''),
                              header = TRUE)
# to keep the primary samples only 
clinical_sample <- read.delim('cBioPortal/pancan_pcawg_2020/data_clinical_sample.txt',
                              header = TRUE, sep = '\t', comment.char = '#', fill = FALSE, quote = '')
primary_sample <- filter(clinical_sample, SAMPLE_TYPE=='Primary') %>%
  inner_join(white_list_sample, by = join_by(SAMPLE_ID==icgc_specimen_id)) %>%
  filter(is.na(PLOIDY) == FALSE)

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

tmb_master <- read.csv(paste(working_dir, '/data/white_list_tmb.csv', sep = ''),
                       header = TRUE)
# TMB of 3 groups, no introns ----
no_introns_smc56_smv <- filter(smc56_smv, Variant_Classification != 'Intron')
cohort_master <- bind_rows(unique(smc56_high_cnv[,c('Tumor_Sample_Barcode','COHORT')]),
                           unique(smc56_low_cnv[,c('Tumor_Sample_Barcode','COHORT')]),
                           unique(no_introns_smc56_smv[,c('Tumor_Sample_Barcode','COHORT')]),
                           unique(smc56_wt))
df <- inner_join(tmb_master, cohort_master, 
                 by=join_by(Tumor_Sample_Barcode), relationship = 'many-to-many') %>%
  mutate(color = with(., ifelse(COHORT %in% c('low CNV', 'high CNV'), 
                                'CNV', COHORT))) %>%
  select(all_of(c('color', 'Tumor_Sample_Barcode', 'TMB', 'abs_mut_number'))) %>%
  unique() %>%
  mutate(color = forcats::fct_relevel(color, 'SmV', 'CNV', 'WT'))
# PLOT 
ggplot(df, aes(x=color, y=log10(TMB))) +
  geom_boxplot() +
  stat_compare_means(method = 'wilcox.test', label = 'p.adj',
                     comparisons = list(c('CNV','WT'), c('SmV','CNV'), c('SmV','WT'))) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_blank()
  ) 
ggsave(paste(working_dir, '/figures/no_introns_3groups_log10tmb_boxplot.eps', sep = ''),
       dpi = 320, width = 6, height = 5, units = 'in')