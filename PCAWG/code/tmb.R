# July 1st, 2025  
# Thi Tran 
# TMB of 3 groups: SmV, CNV, WT

# ENVIRONMENT ####
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
# this file contains the wrong value of TMB since I divided the number of mutations by 30
# so I will recalculate and output it to use for later
tmb_master <- mutate(tmb_master,TMB=abs_mut_number/2800)
write.csv(tmb_master, paste(working_dir,'/data/white_list_tmbv2.csv',sep=''), row.names=F)
# July 1st, 2025 ----
## TMB of 3 groups, INTRON included ----
smc56_wt <- filter(primary_sample, !aliquot_id %in% union(smc56_smv$Tumor_Sample_Barcode,
                                                          union(smc56_low_cnv$Tumor_Sample_Barcode,
                                                                smc56_high_cnv$Tumor_Sample_Barcode))) %>%
  mutate(COHORT = 'WT') %>%
  select(all_of(c('aliquot_id', 'COHORT'))) %>%
  unique()
colnames(smc56_wt) <- c('Tumor_Sample_Barcode', 'COHORT')
cohort_master <- bind_rows(unique(smc56_high_cnv[,c('Tumor_Sample_Barcode','COHORT')]),
                           unique(smc56_low_cnv[,c('Tumor_Sample_Barcode','COHORT')]),
                           unique(smc56_smv[,c('Tumor_Sample_Barcode','COHORT')]), 
                           smc56_wt)
df <- inner_join(tmb_master, cohort_master, 
                by=join_by(Tumor_Sample_Barcode), relationship = 'many-to-many') %>%
  mutate(color = with(., ifelse(COHORT %in% c('low CNV', 'high CNV'), 
                                'CNV', COHORT))) %>%
  select(all_of(c('color', 'Tumor_Sample_Barcode', 'TMB', 'abs_mut_number'))) %>%
  unique() %>%
  mutate(color = forcats::fct_relevel(color, 'SmV', 'CNV', 'WT'))
# PLOT 
ggplot(df, aes(x=color, y=log10(TMB))) +
  geom_violin() + 
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
ggsave(paste(working_dir, '/figures/tmb_3groups_violin.eps', sep = ''),
       dpi = 320, width = 4, height = 4, units = 'in')

ggplot(df, aes(x=color, y=log10(TMB))) +
  geom_boxplot() + 
  stat_compare_means(method = 'wilcox.test', label = 'p.format',
                     comparisons = list(c('CNV','WT'), c('SmV','CNV'), c('SmV','WT'))) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_blank()
  ) 
ggsave(paste(working_dir, '/figures/tmb_3groups_boxplot.eps', sep = ''),
       dpi = 320, width = 4, height = 4, units = 'in')
## TMB of 3 groups, no introns ----
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
  geom_violin() + 
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
ggsave(paste(working_dir, '/figures/no_introns_3groups_log10tmb_violin.eps', sep = ''),
       dpi = 320, width = 6, height = 5, units = 'in')
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
# Sep 5th, 2025 ----
# show everything, plus percentage of tumors with SMC5/6 SmVs
## percentage of tumors with SMC5/6 SmV ####
raw_snv <- read.delim('/storage2/fs1/abby.green/Active/RawData/TCGA_PCAWG/PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf',
                         header = TRUE, sep = '\t', comment.char = '#', fill = FALSE, quote = '')
snv_master <- filter(raw_snv, Tumor_Sample_Barcode %in% white_list_sample$aliquot_id) # only primary samples from white list patients
smc56_smv <- filter(snv_master, Hugo_Symbol %in% smc56_gene) %>%
  filter(Variant_Classification != 'Intron') # only non-intron SMC5/6 SmV
df_plot<- tmb_master %>% 
  filter(Tumor_Sample_Barcode %in% primary_sample$aliquot_id) %>% # only primary samples from white-list patients
  mutate(Tumor_Sample_Barcode = fct_reorder(Tumor_Sample_Barcode, -desc(TMB))) %>%
  mutate(quartile = with(., ifelse(TMB >= 2*3.755446e+00, 'outliers', 
                                   ifelse(TMB <= 7.771429e-01,
                                          'first_quartile','others')))) %>%
  mutate(status = with(., ifelse(Tumor_Sample_Barcode %in% smc56_smv$Tumor_Sample_Barcode,
                                 'SMC5/6 non-intronic SmV', 'No SMC5/6 non-intronic SmV')))
df_plot %>% group_by(quartile, status) %>%
  summarise(n=n())
# low_quartile = 100*11/(11+433) = 2.477%
# others = 100*124/(124+1068) = 10.403%
# outliers = 100*47/(47+91) = 34.058%
## highlight the quartile -----
df_plot %>%
  ggplot(aes(x=Tumor_Sample_Barcode,y=TMB,color=quartile)) +
  geom_point() +
  geom_vline(xintercept = c('fc8130e0-ad66-b82e-e040-11ac0d485e0e','48624a82-c623-11e3-bf01-24c6515278c0'),
             linetype = 'longdash') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = 'top') +
  scale_color_manual(values = c('#5773CCFF','#5773CCFF','#FFB900FF'))
ggsave(paste(working_dir, '/figures/all_pass_tmb.eps', sep = ''), dpi = 320)
## highlight the SMC5/6 tumors -----
df_plot %>%
  ggplot(aes(x=Tumor_Sample_Barcode,y=TMB,color=status)) +
  geom_point() +
  geom_vline(xintercept = c('fc8130e0-ad66-b82e-e040-11ac0d485e0e','48624a82-c623-11e3-bf01-24c6515278c0'),
             linetype = 'longdash') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = 'top') +
  scale_color_manual(values = c('black', 'red'))
ggsave(paste(working_dir, '/figures/highlight_all_pass_tmb.eps', sep = ''), dpi = 320)
# Oct 29th, 2025 ----
# intronic vs exonic SMC5/6 SmVs 
smc56_intron_mut <- filter(snv_master, Hugo_Symbol %in% smc56_gene & Variant_Classification =='Intron') %>%
  mutate(COHORT='intron SmV')
cohort_master2 <- rbind(unique(smc56_intron_mut[,c('Tumor_Sample_Barcode','COHORT')]),
                        cohort_master)
df_plot <- filter(tmb_master, Tumor_Sample_Barcode %in% primary_sample$aliquot_id) %>%
  inner_join(cohort_master2, by=join_by(Tumor_Sample_Barcode), relationship='one-to-many') %>%
  mutate(color=with(., ifelse(COHORT %in% c('low CNV','high CNV'), 'CNV', COHORT))) %>%
  select(all_of(c('color','Tumor_Sample_Barcode','TMB','abs_mut_number'))) %>%
  unique() %>%
  mutate(color=fct_relevel(color,'intron SmV','SmV','CNV','WT'))
ggplot(df_plot, aes(x=color,y=log10(TMB))) +
  geom_boxplot() +
  stat_compare_means(method = 'wilcox.test', label = 'p.format',
                     comparisons = list(c('intron SmV','SmV'),c('SmV','CNV'),
                                        c('CNV','WT'),c('intron SmV','CNV'),
                                        c('SmV','WT'),c('intron SmV','WT'))) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_blank()
  ) 
ggsave(paste(working_dir,'/figures/log10TMB_intron_smc56_smv_boxplot.eps', sep=''),
       dpi = 320, width = 5, height = 6.5, units = 'in')

ggplot(df_plot, aes(x=color,y=log10(TMB))) +
  geom_violin() +
  stat_compare_means(method = 'wilcox.test', label = 'p.format',
                     comparisons = list(c('intron SmV','SmV'),c('SmV','CNV'),
                                        c('CNV','WT'),c('intron SmV','CNV'),
                                        c('SmV','WT'),c('intron SmV','WT'))) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_blank()
  ) 
ggsave(paste(working_dir,'/figures/log10TMB_intron_smc56_smv_violin.eps', sep=''),
       dpi = 320, width = 5, height = 6.5, units = 'in')
