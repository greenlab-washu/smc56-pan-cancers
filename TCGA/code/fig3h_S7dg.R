# Fig3, FigS7 - d, g
# ENVIRONMENT ----
library('tidyverse')
library('ggplot2')
library('survival')
library('ggpubr')
library('survminer')
input_dir <- 'TCGA'
output_dir <- 'TCGA'
smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
                'NSMCE4A', 'NSMCE4B', 'EID3',
                'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')
theme_general <- theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

# TCGA-UCEC and TCGA-UCS ----

## INPUT ----
snv_master <- read.csv(paste(input_dir, '/data/gdc_snv.csv', sep = ''),
                       header = T)
tmb_master <- read.csv(paste(input_dir, '/data/gdc_tmb.csv', sep = ''),
                       header = T)
cohort_sample <- 
  list.files(path = paste(input_dir,'/data',sep=''),
             full.names = TRUE, pattern = 'cohort_for_fpkm.csv$') %>%
  purrr::map(read.csv, header = TRUE) 
colnames(cohort_sample[[2]]) <- colnames(cohort_sample[[4]])
cohort_uterus <- bind_rows(cohort_sample) %>%
  filter(project %in% c('TCGA-UCEC','TCGA-UCS')) %>% 
  mutate(Patient_Barcode=str_extract(Tumor_Sample_Barcode,'^.{12}'))

ucec_survival <- read.delim(paste(input_dir, '/cBioPortal/pan_can_atlas/ucec_tcga_pan_can_atlas_2018/data_clinical_patient.txt',
                                  sep=''),
                            sep = '\t', header = T, comment.char = '#') %>%
  mutate_at(c('DAYS_LAST_FOLLOWUP','DAYS_TO_BIRTH','AJCC_STAGING_EDITION'), as.integer)
ucs_survival <- read.delim(paste(input_dir, '/cBioPortal/pan_can_atlas/ucs_tcga_pan_can_atlas_2018/data_clinical_patient.txt',
                                  sep=''),
                            sep = '\t', header = T, comment.char = '#') %>%
  mutate_at(c('DAYS_LAST_FOLLOWUP','DAYS_TO_BIRTH','AJCC_STAGING_EDITION'), as.integer)
cbioportal_clinical <- bind_rows(ucec_survival,ucs_survival) %>%
  mutate(OS_STATUS_DIGIT=as.numeric(str_extract(OS_STATUS,'^.{1}'))) %>%
  filter(!is.na(OS_MONTHS))

ucec_ploidy <- read.delim(paste(input_dir, '/cBioPortal/pan_can_atlas/ucec_tcga_pan_can_atlas_2018/data_clinical_sample.txt',
                                  sep=''),
                            sep = '\t', header = T, comment.char = '#')
ucs_ploidy <- read.delim(paste(input_dir, '/cBioPortal/pan_can_atlas/ucs_tcga_pan_can_atlas_2018/data_clinical_sample.txt',
                                 sep=''),
                           sep = '\t', header = T, comment.char = '#') 
cbioportal_ploidy <- bind_rows(ucec_ploidy,ucs_ploidy) 

## 3h - 3 groups ----
df <- left_join(cohort_uterus, tmb_master, 
                by=join_by(Tumor_Sample_Barcode), relationship='many-to-one') %>%
  inner_join(cbioportal_clinical,by=join_by(Patient_Barcode==PATIENT_ID)) %>%
  mutate(status=with(.,ifelse(COHORT=='wt','WT',
                              ifelse(COHORT %in% c('low_cnv','high_cnv'),'CNV','SmV')))) %>%
  mutate(status=fct_relevel(status, 'SmV','CNV','WT')) %>%
  mutate(COHORT=fct_relevel(COHORT,'syn','del','nonsyn_nondel',
                            'low_cnv','high_cnv','wt'))
# count the number of tumors in each group
summarise(unique(df[,c('Tumor_Sample_Barcode','status')]),.by=status,n=n())

select(df,all_of(c('Tumor_Sample_Barcode','TMB','status'))) %>%
  unique() %>%
  ggplot(aes(x=status, y=log10(TMB), color=status)) +
  geom_boxplot() +
  theme_general + 
  theme(legend.position = 'none') +
  scale_color_manual(values=c('#A52A2A','#702963','#28446f')) +
  stat_compare_means(comparisons = list(c('SmV','CNV'),c('CNV','WT'),c('SmV','WT')),
                     method = 'wilcox', label = 'p.adjust')
ggsave(paste(output_dir,'/TCGA/figures/uterus_log10TMB_3groups.eps',sep=''),
       dpi=320, width=3, height=3, units='in')

## S7d - all subgroups ----
# count the number of tumors in each group 
summarise(unique(df[,c('Tumor_Sample_Barcode','COHORT')]),.by=COHORT,n=n())
df %>% select(all_of(c('Tumor_Sample_Barcode','TMB','COHORT'))) %>%
  unique() %>%
  ggplot(aes(x=COHORT, y=log10(TMB), color=COHORT)) +
  geom_boxplot() +
  theme_general + 
  theme(legend.position = 'none') +
  scale_color_manual(values=c('#FAA0A0', '#A52A2A', '#80B19B', '#C6B3DD', '#702963', '#28446f')) +
  stat_compare_means(comparisons = list(c('syn','nonsyn_nondel'), c('nonsyn_nondel','del'), c('syn','del'),
                                        c('syn','wt'), c('nonsyn_nondel','wt'), c('del','wt'),
                                        c('low_cnv','wt'),c('high_cnv','wt')),
                     method = 'wilcox', label = 'p.adjust')
ggsave(paste(output_dir,'/TCGA/figures/uterus_log10TMB_6groups.eps',sep=''),
       dpi=320, width=4.5, height=3, units='in')

## ANEUPLOIDY ----
### S7g - 3 main groups -----
df_ploidy <- inner_join(df,cbioportal_ploidy,
                        by=join_by(Patient_Barcode==PATIENT_ID),relationship = 'many-to-one') %>%
  mutate(COHORT=fct_relevel(COHORT,'syn','del','nonsyn_nondel',
                            'low_cnv','high_cnv','wt')) %>%
  mutate(status=fct_relevel(status,'SmV','CNV','WT'))
unique(df_ploidy[,c('Tumor_Sample_Barcode','status')]) %>%
  summarise(.by=status,n=n())
select(df_ploidy,all_of(c('Tumor_Sample_Barcode','status','ANEUPLOIDY_SCORE'))) %>%
  unique() %>%
  ggplot(aes(x=status,y=ANEUPLOIDY_SCORE,color=status)) +
  geom_boxplot() +
  theme_general +
  theme(legend.position = 'none') +
  scale_color_manual(values=c('#A52A2A','#702963','#28446f')) +
  stat_compare_means(comparisons = list(c('SmV','CNV'),c('CNV','WT'),c('SmV','WT')),
                     method = 'wilcox',label='p.adjust')
ggsave(paste(output_dir,'/TCGA/figures/uterus_aneuploidy_3groups.eps',sep=''),
       dpi=320, width=3, height=3, units='in')