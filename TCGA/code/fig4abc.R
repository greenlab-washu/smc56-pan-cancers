# ENVIRONMENT ----
install.packages('survminer')
library('survminer')
require('survival')
library('ggplot2')
library('tidyverse')

inout_dir <- 'TCGA'
output_dir <- 'TCGA'

gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
          'NSMCE4A', 'NSMCE4B', 'EID3',
          'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')
group_color <- c('#FAA0A0', '#A52A2A', '#80B19B', '#C6B3DD', '#702963', '#28446f')
groups <- c('syn', 'del', 'nonsyn_nondel', 'low_cnv', 'high_cnv', 'wt')
theme_general <- theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())
# INPUT ---- 
## survival files ----
patient_survival_files <- list.files(path = dir('cBioPortal/pan_can_atlas',
                                        full.names = TRUE),
                             full.names = TRUE, pattern = 'data_clinical_patient.txt')
patient_tcga_survival <- purrr::map(patient_survival_files, read.delim,
                                    sep = '\t', header = TRUE, comment.char = '#', fill = FALSE) %>%
  purrr::map(mutate_at, c('DAYS_LAST_FOLLOWUP', 'DAYS_TO_BIRTH', 'AJCC_STAGING_EDITION'), as.integer)
tcga_survival <- bind_rows(patient_tcga_survival)
## cohort ----
cohort_sample <- 
  list.files(path = input_dir, full.names = TRUE, pattern = 'cohort_for_fpkm.csv$') %>%
  purrr::map(read.csv, header = TRUE) 
colnames(cohort_sample[[2]]) <- colnames(cohort_sample[[4]])
cohort_sample <- bind_rows(cohort_sample)

cohort_sample$Patient_Barcode <- stringr::str_extract(
  cohort_sample$Tumor_Sample_Barcode, '^.{12}'
)
# ANALYSIS ----
df <- 
  merge(tcga_survival, cohort_sample, by.x = 'PATIENT_ID', by.y = 'Patient_Barcode') %>%
  mutate(OS_STATUS = as.integer(stringr::str_extract(OS_STATUS, '[:digit:]'))) %>%
  mutate(COHORT = forcats::fct_relevel(COHORT, 'syn', 'del', 'nonsyn_nondel', 
                                       'low_cnv', 'high_cnv', 'wt')) %>%
  select(any_of(c('COHORT', 'PATIENT_ID', 'OS_MONTHS', 'OS_STATUS'))) %>%
  unique() %>%
  mutate(OS_STATUS=as.numeric(as.integer(stringr::str_extract(OS_STATUS, '[:digit:]'))))

## 4a - altered vs WT ----
df_altered <- df %>%
  mutate(group=with(.,ifelse(COHORT=='wt','altered','WT'))) %>%
  select(all_of(c('OS_MONTHS','OS_STATUS','group'))) %>%
  unique() %>%
  mutate(group=fct_relevel(group,'altered','WT'))
fit <- surv_fit(Surv(OS_MONTHS,OS_STATUS) ~ group, df_altered)
df_altered %>%
  ggsurvplot(fit, data = ., pval = T,
             ggtheme = theme_general, legend = 'right',
             legend.title = '', xlab = 'Time (month)')
ggsave(paste(output_dir, '/figures/altered_survival.eps', sep = ''),
       dpi = 320, width = 5, height = 3)

#the pairwise comparison 
pairwise_survdiff(Surv(OS_MONTHS, OS_STATUS) ~ group, data = df_altered,
                  p.adjust.method = 'none')

## 4b - 3 groups ----
df_3groups <- df %>%
  mutate(group=ifelse(with(., COHORT=='wt','wt',
                           ifelse(COHORT %in% c('low_cnv','high_cnv'),
                                  'CNA','SmV')))) %>%
  select(all_of(c('group','OS_MONTHS','OS_STATUS'))) %>%
  unique() %>%
  mutate(group=fct_relevel(group,'SmV','CNV','wt'))
fit <- surv_fit(Surv(OS_MONTHS, OS_STATUS) ~ group, data = df_3groups)
df_3groups %>% 
  ggsurvplot(fit, data = ., pval = TRUE,
             palette = c('#A52A2A', '#702963', '#28446f'), 
             ggtheme = theme_general,
             legend = 'right', legend.title = '',
             xlab = 'Time (month)')
ggsave(paste(output_dir, '/figures/3groups_survival.eps', sep = ''),
       dpi = 320, width = 5, height = 3)
#the pairwise comparison 
pairwise_survdiff(Surv(OS_MONTHS, OS_STATUS) ~ group, data = df_3groups,
                  p.adjust.method = 'none')


## 4c - different types of SmVs vs WT ----
fit <- filter(df, COHORT %in% c('syn', 'del', 'nonsyn_nondel', 'wt')) %>%
  surv_fit(Surv(OS_MONTHS, OS_STATUS) ~ COHORT, data = .)
df %>% filter(COHORT %in% c('syn', 'del', 'nonsyn_nondel', 'wt')) %>%
  mutate(OS_STATUS = as.integer(stringr::str_extract(OS_STATUS, '[:digit:]'))) %>%
  ggsurvplot(fit, data = ., pval = TRUE,
             palette = c('#FAA0A0', '#A52A2A', '#80B19B', '#28446f'), 
             ggtheme = theme_general,
             legend = 'right', legend.title = '',
             xlab = 'Time (month)')
ggsave(paste(output_dir, '/figures/separated_snv_abs_cnv_survival.eps', sep = ''),
       dpi = 320, width = 5, height = 3)
#the pairwise comparison 
filter(df, COHORT %in% c('syn','del','nonsyn_nondel','wt')) %>%
  pairwise_survdiff(Surv(OS_MONTHS, OS_STATUS) ~ COHORT, data = .,
                  p.adjust.method = 'none')
