# May 6th, 2025
# Thi Tran 
# Survival analysis for Hartwig patients 
# (metastsis colorectal samples)

# ENVIRONMENT ----
library('tidyverse')
library('lubridate')
library('survival')
library('ggplot2')
#install.packages(c('ggsurvfit', 'gtsummary'))
library('ggsurvfit')
library('gtsummary')
library('ggpubr',
        lib.loc = '/storage2/fs1/abby.green/Active/R_libraries/Bunny-Wunnies Freak Out.Bunny-Wunnies Freak Out')
library('survminer',
        lib.loc = '/storage2/fs1/abby.green/Active/R_libraries/Bunny-Wunnies Freak Out.Bunny-Wunnies Freak Out')
theme_general <- theme_classic() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

input_dir <- '/storage2/fs1/abby.green/Active/Hartwig'
output_dir <- '/storage2/fs1/abby.green/Active/Users/thi/SMC56_de_novo/2025/Hartwig'

clinical <- read_tsv(paste(input_dir, '/metadata.tsv', sep = ''), na = c('null'))
post_biopsy <- read_tsv(paste(input_dir, '/post_biopsy_drugs.tsv', sep = ''), na = c('null'))
pre_biopsy <- read_tsv(paste(input_dir, '/pre_biopsy_drugs.tsv', sep = ''), na = c('null'))
# August 1st, 2025 - omit the patients with missing treatment response information ####
## INPUT ----
treatment <- read_tsv(paste(input_dir, '/treatment_responses.tsv', sep = ''), na = c('null')) %>%
  group_by(patientId) %>% 
  filter(any(!is.na(response))) %>% # omit patients with no response record
  ungroup()
cohort_files <- list.files(paste(output_dir, '/data', sep = ''),
                           pattern = '_smc56.csv$', full.names = TRUE)
cohort_indv <- purrr::map(cohort_files, read.csv, header = TRUE)
# master survival timeline
master_survival <- inner_join(clinical[,c('sampleId','hmfPatientId','gender','deathDate','biopsyDate',
                             'hasSystemicPreTreatment','hasRadiotherapyPreTreatment')],
                 pre_biopsy[,c('sampleId','startDate')], # use the earliest start date of this data table as the beginning of it all
                 by=join_by(sampleId), relationship='one-to-many') %>%
  inner_join(treatment[,c('sampleId','startDate','responseDate')], # for the end date of patients with missing death date
             by=join_by(sampleId), relationship='many-to-many') %>%
  mutate(status = with(., ifelse(is.na(deathDate),0,1)))
colnames(master_survival) <- c('sampleId','hmfPatientId','gender','deathDate','biopsyDate',
                  'hasSystemicPreTreatment','hasRadiotherapyPreTreatment',
                  'preBiopsyStartDate','postBiopsyStartDate','postBiopsyResponseDate','status')
condense_survival <- group_by(master_survival, sampleId) %>%
  mutate(lastResponseDate = max(postBiopsyResponseDate)) %>%
  ungroup() %>%
  mutate(survival_time = with(., ifelse(is.na(preBiopsyStartDate), 
                                        ifelse(is.na(deathDate),
                                                     as.duration(postBiopsyStartDate %--% lastResponseDate)/dmonths(1),
                                                     as.duration(postBiopsyStartDate %--% deathDate)/dmonths(1)),
                                        ifelse(is.na(deathDate),
                                               as.duration(preBiopsyStartDate %--% lastResponseDate)/dmonths(1),
                                               as.duration(preBiopsyStartDate %--% deathDate)/dmonths(1))))) %>%
  group_by(sampleId) %>%
  mutate(time = max(survival_time)) %>% # there could be many pre biopsy start date, which leads to different survival time for the same patient
  ungroup() %>%
  select(all_of(c('sampleId','gender','time','status',
                  'hasSystemicPreTreatment','hasRadiotherapyPreTreatment'))) %>%
  unique()
## GRAPHS ----
# DEL vs WT
df_del <- filter(condense_survival, sampleId %in% cohort_indv[[3]]$Tumor_Sample_Barcode) %>%
  mutate(COHORT = 'Del SnV')
df_wt <- filter(condense_survival, sampleId %in% cohort_indv[[6]]$Tumor_Sample_Barcode) %>%
  mutate(COHORT = 'WT')
df <- rbind(df_del, df_wt) 
fit <- survfit(Surv(time, status) ~ COHORT, data = df)
ggsurvplot(fit, data = df, pval = TRUE, 
           palette = c('#A52A2A','#28446f'), ggtheme = theme_general,
           legend = 'right', legend.title = '',
           xlab = 'Time (month)') 
ggsave(paste(output_dir, '/figures/survival_del_vs_wt.eps', sep = ''),
       dpi = 320, width = 5, height = 3, units = 'in')
coxph(Surv(time, status) ~ COHORT, data = df) 
# SmV vs Everything else
df_smv <- filter(condense_survival, sampleId %in% union(cohort_indv[[3]]$Tumor_Sample_Barcode,
                                                        union(cohort_indv[[4]]$Tumor_Sample_Barcode, cohort_indv[[5]]$Tumor_Sample_Barcode))) %>%
  mutate(COHORT = 'SmV')
df_others <- filter(condense_survival, sampleId %in% union(cohort_indv[[6]]$Tumor_Sample_Barcode,
                                                           union(cohort_indv[[1]]$Tumor_Sample_Barcode, cohort_indv[[2]]$Tumor_Sample_Barcode))) %>%
  mutate(COHORT = 'others')
df <- rbind(df_smv, df_others)
fit <- survfit(Surv(time, status) ~ COHORT, data = df)
ggsurvplot(fit, data = df, pval = TRUE, 
           palette = c('#A52A2A','#28446f'), ggtheme = theme_general,
           legend = 'right', legend.title = '',
           xlab = 'Time (month)') 
ggsave(paste(output_dir, '/figures/survival_smv_others.eps', sep = ''),
       dpi = 320, width = 5, height = 3, units = 'in')
coxph(Surv(time, status) ~ COHORT, data = df) 
# coding, non-coding, CNV, and WT
smc56_snv <- read.csv(paste(output_dir, '/data/smc56_snv.csv', sep = ''),
                                 header = TRUE, na.strings = c('NA', ''))
no_introns_smc56_snv <- filter(smc56_snv, Variant_Classification != 'Intron')
coding_smc56 <-filter(no_introns_smc56_snv, !is.na(HGVSp_Short))
non_coding_smc56 <- filter(no_introns_smc56_snv, is.na(HGVSp_Short) & !Tumor_Sample_Barcode %in% coding_smc56$Tumor_Sample_Barcode )
df_coding <- filter(condense_survival, sampleId %in% coding_smc56$Tumor_Sample_Barcode) %>%
  mutate(COHORT = 'Coding SnV')
df_non_coding <- filter(condense_survival, sampleId %in% non_coding_smc56$Tumor_Sample_Barcode) %>%
  mutate(COHORT = 'Non-coding SnV')
df_cnv <- filter(condense_survival, sampleId %in% union(cohort_indv[[1]]$Tumor_Sample_Barcode, cohort_indv[[2]]$Tumor_Sample_Barcode)) %>%
  mutate(COHORT = 'CNV')
df <- rbind(df_coding, df_non_coding) %>%
  rbind(df_cnv) %>%
  rbind(df_wt)
fit <- survfit(Surv(time, status) ~ COHORT, data = df)
ggsurvplot(fit, data = df, pval = TRUE, 
           palette = c('#702963','#A52A2A','#80B19B','#28446f'), 
           ggtheme = theme_general,
           legend = 'right', legend.title = '',
           xlab = 'Time (month)')
ggsave(paste(output_dir, '/figures/coding_non_coding_survival.eps', sep = ''),
       dpi = 320, width = 5, height = 3, units = 'in')
pairwise_survdiff(Surv(time, status) ~ COHORT, data = df)
# August 13th - include all patients - deal with missing information later ####
## INPUT ----
treatment <- read_tsv(paste(input_dir, '/treatment_responses.tsv', sep = ''), na = c('null')) 
cohort_files <- list.files(paste(output_dir, '/data', sep = ''),
                           pattern = 'patients_.+_smc56.csv$', full.names = TRUE)
cohort_indv <- purrr::map(cohort_files, read.csv, header = TRUE)
# SmV
del_indv <- rbind(cohort_indv[[1]], cohort_indv[[5]]) %>%
  rbind(cohort_indv[[10]])
nsnd_indv <- rbind(cohort_indv[[2]], cohort_indv[[6]]) %>%
  rbind(cohort_indv[[11]])
syn_indv <- cohort_indv[[3]]
# CNV
low_cnv_indv <- cohort_indv[[9]]
high_cnv_indv <- cohort_indv[[8]]
# WT
wt_indv <- rbind(cohort_indv[[4]], cohort_indv[[7]]) %>%
  rbind(cohort_indv[[12]])
# master survival timeline
master_survival <- inner_join(clinical[,c('sampleId','hmfPatientId','gender','deathDate','biopsyDate',
                                          'hasSystemicPreTreatment','hasRadiotherapyPreTreatment')],
                              pre_biopsy[,c('sampleId','startDate')], # use the earliest start date of this data table as the beginning of it all
                              by=join_by(sampleId), relationship='one-to-many') %>%
  inner_join(treatment[,c('sampleId','startDate','responseDate')], # for the end date of patients with missing death date
             by=join_by(sampleId), relationship='many-to-many') %>%
  mutate(status = with(., ifelse(is.na(deathDate),0,1)))
colnames(master_survival) <- c('sampleId','hmfPatientId','gender','deathDate','biopsyDate',
                               'hasSystemicPreTreatment','hasRadiotherapyPreTreatment',
                               'preBiopsyStartDate','postBiopsyStartDate','postBiopsyResponseDate','status')
condense_survival <- group_by(master_survival, sampleId) %>%
  mutate(lastResponseDate = max(postBiopsyResponseDate)) %>%
  ungroup() %>%
  mutate(survival_time = with(., ifelse(is.na(preBiopsyStartDate), 
                                        ifelse(is.na(deathDate),
                                               as.duration(postBiopsyStartDate %--% lastResponseDate)/dmonths(1),
                                               as.duration(postBiopsyStartDate %--% deathDate)/dmonths(1)),
                                        ifelse(is.na(deathDate),
                                               as.duration(preBiopsyStartDate %--% lastResponseDate)/dmonths(1),
                                               as.duration(preBiopsyStartDate %--% deathDate)/dmonths(1))))) %>%
  group_by(sampleId) %>%
  mutate(time = max(survival_time)) %>% # there could be many pre biopsy start date, which leads to different survival time for the same patient
  ungroup() %>%
  select(all_of(c('sampleId','gender','time','status',
                  'hasSystemicPreTreatment','hasRadiotherapyPreTreatment'))) %>%
  unique()
# output this file so I can quickly filter out sample with missing information 
# in downstream analysis 
write.csv(condense_survival, paste(output_dir, '/data/survival.csv', sep = ''), row.names = FALSE)
## GRAPHS ----
# DEL vs WT - there is no overlapping samples between these two groups
df_del <- filter(condense_survival, !is.na(time)) %>%
  filter(sampleId %in% del_indv$Tumor_Sample_Barcode) %>%
  mutate(COHORT = 'Del SnV')
df_wt <- filter(condense_survival, !is.na(time)) %>%
  filter(sampleId %in% wt_indv$Tumor_Sample_Barcode) %>%
  mutate(COHORT = 'WT')
df <- rbind(df_del, df_wt) 
fit <- survfit(Surv(time, status) ~ COHORT, data = df)
ggsurvplot(fit, data = df, pval = TRUE, 
           palette = c('#A52A2A','#28446f'), ggtheme = theme_general,
           legend = 'right', legend.title = '',
           xlab = 'Time (month)') 
ggsave(paste(output_dir, '/figures/all_patients_survival_del_vs_wt.eps', sep = ''),
       dpi = 320, width = 5, height = 3, units = 'in')
coxph(Surv(time, status) ~ COHORT, data = df) 
# SmV vs everything else - there are overlapping tumors in these 2 groups
df_smv <- filter(condense_survival, !is.na(time) &
                   sampleId %in% union(del_indv$Tumor_Sample_Barcode,
                                       union(nsnd_indv$Tumor_Sample_Barcode, syn_indv$Tumor_Sample_Barcode))) %>%
  mutate(COHORT = 'SmV')
df_others <- filter(condense_survival, !is.na(time) & 
                      sampleId %in% union(wt_indv$Tumor_Sample_Barcode,
                                          union(low_cnv_indv$Tumor_Sample_Barcode, high_cnv_indv$Tumor_Sample_Barcode))) %>%
  mutate(COHORT = 'others')
df <- rbind(df_smv, df_others)
fit <- survfit(Surv(time, status) ~ COHORT, data = df)
ggsurvplot(fit, data = df, pval = TRUE, 
           palette = c('#A52A2A','#28446f'), ggtheme = theme_general,
           legend = 'right', legend.title = '',
           xlab = 'Time (month)') 
ggsave(paste(output_dir, '/figures/all_patients_survival_smv_others.eps', sep = ''),
       dpi = 320, width = 5, height = 3, units = 'in')
coxph(Surv(time, status) ~ COHORT, data = df) 
# coding, non-coding, CNV, and WT - there are overlapping tumors between coding, non-coding, and CNV
smc56_snv_files <- list.files(paste(output_dir, '/data', sep = ''), 
                              pattern = 'patients_smc56_snv.csv$', full.names = TRUE) 
smc56_snv_indv <- purrr::map(smc56_snv_files, read.csv, header = TRUE, na.strings = c('NA', '')) 
smc56_snv <- rbind(smc56_snv_indv[[1]], smc56_snv_indv[[2]], smc56_snv_indv[[3]])
rm(smc56_snv_files, smc56_snv_indv)
no_introns_smc56_snv <- filter(smc56_snv, Variant_Classification != 'Intron')
coding_smc56 <-filter(no_introns_smc56_snv, !is.na(HGVSp_Short))
non_coding_smc56 <- filter(no_introns_smc56_snv, is.na(HGVSp_Short) & !Tumor_Sample_Barcode %in% coding_smc56$Tumor_Sample_Barcode )
df_coding <- filter(condense_survival, !is.na(time) & 
                      sampleId %in% coding_smc56$Tumor_Sample_Barcode) %>%
  mutate(COHORT = 'Coding SnV')
df_non_coding <- filter(condense_survival, !is.na(time) &
                          sampleId %in% non_coding_smc56$Tumor_Sample_Barcode) %>%
  mutate(COHORT = 'Non-coding SnV')
df_cnv <- filter(condense_survival, !is.na(time) &
                   sampleId %in% union(low_cnv_indv$Tumor_Sample_Barcode, high_cnv_indv$Tumor_Sample_Barcode)) %>%
  mutate(COHORT = 'CNV')
df <- rbind(df_coding, df_non_coding) %>%
  rbind(df_cnv) %>%
  rbind(df_wt)
fit <- survfit(Surv(time, status) ~ COHORT, data = df)
ggsurvplot(fit, data = df, pval = TRUE, 
           palette = c('#702963','#A52A2A','#80B19B','#28446f'), 
           ggtheme = theme_general,
           legend = 'right', legend.title = '',
           xlab = 'Time (month)')
ggsave(paste(output_dir, '/figures/all_patients_coding_non_coding_survival.eps', sep = ''),
       dpi = 320, width = 5, height = 3, units = 'in')
pairwise_survdiff(Surv(time, status) ~ COHORT, data = df)
# September 2nd - biallelic patients ----
smc56_biallelic <- read.csv(paste(output_dir, '/data/smc56_biallelic.csv', sep = ''), header = TRUE)
# analyze only those patients with survival information 
condense_survival <- read.csv(paste(output_dir, '/data/survival.csv', sep = ''), header = TRUE)
## INPUT ----
treatment <- read_tsv(paste(input_dir, '/treatment_responses.tsv', sep = ''), na = c('null')) 
cohort_files <- list.files(paste(output_dir, '/data', sep = ''),
                           pattern = 'patients_.+_smc56.csv$', full.names = TRUE)
cohort_indv <- purrr::map(cohort_files, read.csv, header = TRUE)
# SmV
del_indv <- rbind(cohort_indv[[1]], cohort_indv[[5]]) %>%
  rbind(cohort_indv[[10]])
nsnd_indv <- rbind(cohort_indv[[2]], cohort_indv[[6]]) %>%
  rbind(cohort_indv[[11]])
syn_indv <- cohort_indv[[3]]
# CNV
low_cnv_indv <- cohort_indv[[9]]
high_cnv_indv <- cohort_indv[[8]]
# WT
wt_indv <- rbind(cohort_indv[[4]], cohort_indv[[7]]) %>%
  rbind(cohort_indv[[12]])

## no intron biallelic SmV, non-biallelic SmV, vs others ----
no_intron_smc56_biallelic <- filter(smc56_biallelic, Variant_Classification != 'Intron')
df <- filter(condense_survival, !is.na(time)) %>%
  mutate(group = with(., ifelse(sampleId %in% no_intron_smc56_biallelic$Tumor_Sample_Barcode,
                                'biallelic SMC5/6', ifelse(sampleId %in% union(del_indv$Tumor_Sample_Barcode,
                                                                               union(nsnd_indv$Tumor_Sample_Barcode,
                                                                                     syn_indv$Tumor_Sample_Barcode)),
                                                           'non-biallelic SMC5/6', 'others'))))

fit <- survfit(Surv(time, status) ~ group, data = df)
ggsurvplot(fit, data = df, pval = TRUE, 
           palette = c('#B84E5AFF','#D8ADC1FF','#28446f'), ggtheme = theme_general,
           legend = 'right', legend.title = '',
           xlab = 'Time (month)') 
ggsave(paste(output_dir, '/figures/all_patients_survival_no_introns_biallelic_3groups.eps', sep = ''),
       dpi = 320, width = 5, height = 3, units = 'in')
pairwise_survdiff(Surv(time, status) ~ group, data = df) # no p-value < 0.05
## biallelic del SMC56 SmV ----
del_smc56_biallelic <- filter(smc56_biallelic, IMPACT == 'HIGH' | 
                                stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') <= 0.05)
df <- filter(condense_survival, !is.na(time)) %>%
  mutate(group = with(., ifelse(sampleId %in% del_smc56_biallelic$Tumor_Sample_Barcode,
                                'del biallelic SmV', ifelse(sampleId %in% setdiff(no_intron_smc56_biallelic$Tumor_Sample_Barcode, del_smc56_biallelic$Tumor_Sample_Barcode),
                                                            'non-del biallelic SmV', ifelse(sampleId %in% union(del_indv$Tumor_Sample_Barcode,
                                                                                                                union(nsnd_indv$Tumor_Sample_Barcode, syn_indv$Tumor_Sample_Barcode)),
                                                                                            'non-biallelic SmV', 'others')))))
# there are no patients with deleterious biallelic SMC5/6 SmV and survival information 
# so abort this mission here 