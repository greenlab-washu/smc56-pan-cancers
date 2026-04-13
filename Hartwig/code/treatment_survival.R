# August 13th, 2025 
# Thi Tran 
# SMC5/6 dysfunction, types of treatments, and patients survival
# August 13th - I think the post biopsy treatment is more appropriate in this context

# ENVIRONMENT ####
library('tidyverse')
library('lubridate')
library('survival')
library('ggplot2')
library('ggsurvfit')
library('gtsummary')
library('ggVennDiagram')
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

smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
                'NSMCE4A', 'NSMCE4B', 'EID3',
                'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')

# August 13th - include all patients, then omit them once they have NA as calculated survival time 
# INPUT ####
clinical <- read_tsv(paste(input_dir, '/metadata.tsv', sep = ''), na = c('null', ''))
post_biopsy <- read_tsv(paste(input_dir, '/post_biopsy_drugs.tsv', sep = ''), na = c('null', ''))
pre_biopsy <- read_tsv(paste(input_dir, '/pre_biopsy_drugs.tsv', sep = ''), na = c('null', ''))
treatment <- read_tsv(paste(input_dir, '/treatment_responses.tsv', sep = ''), na = c('null', '')) 
# load the new cohort files generated on Oct 28th
cohort_master <- list.files(paste(output_dir, '/data', sep = ''),
                            full.names = T, pattern = 'cohort.csv$') %>%
  map(read.csv, header = T) %>%
  bind_rows()
condense_survival <- read.csv(paste(output_dir, '/data/survival.csv', sep = ''),
                              header = TRUE, na.strings = c('NA', ''))
tmb_master <- read.csv(paste(output_dir, '/data/pass_tmb_master.csv', sep = ''),
                       header = TRUE)
pass_smv_master <- list.files(paste(output_dir, '/data', sep = ''), 
                              full.names = T, pattern = 'patients_pass_only') %>% 
  map(read.csv, header=T) %>% 
  bind_rows()
# August 18th - based on post-biopsy treatment
# OVERAL survival time ----
## treatment only ----
df <- inner_join(condense_survival, post_biopsy, 
                 by = join_by(sampleId), relationship = 'one-to-many') %>%
  select(all_of(c('patientId','time','status','type'))) %>%
  unique() %>% 
  filter(!is.na(time)) 
fit <- survfit(Surv(time, status) ~ type, data = df)
pairwise_survdiff(Surv(time, status) ~ type, data = df) 
ggsurvplot(fit, data = df, pval = FALSE, ggtheme = theme_general, 
           palette = c('#C3A016FF','#C3D878FF','#58A787FF','#8EBACDFF','#246893FF','#163274FF','#0C1F4BFF'),
           legend = 'bottom', legend.title = '', 
           xlab = 'Overall survival time (month)', xlim=c(0,150)) 
ggsave(paste(output_dir, '/figures/survival_treatment_overal_time.eps', sep = ''),
       dpi = 320, width = 5.5, height = 4, units = 'in')
## treatment only, with SMC5/6 SmVs ----
smc56_smv_tumors <- filter(cohort_master, Hugo_Symbol %in% smc56_gene &
                             alteration_type %in% c('del_mutation','nsnd_mutation','syn_mutation'))
df <- inner_join(condense_survival, post_biopsy, 
                 by = join_by(sampleId), relationship = 'one-to-many') %>%
  filter(!is.na(time)) %>%
  unique() %>% 
  mutate(group = with(., ifelse(sampleId %in% smc56_smv_tumors$Tumor_Sample_Barcode,
                                'SMC56 SmV','others'))) %>%
  select(all_of(c('patientId','time','status','type','group'))) %>%
  unique() %>% 
  filter(!is.na(type)) %>%
  mutate(group=fct_relevel(group, 'SMC56 SmV', 'others')) %>%
  mutate(type = fct_relevel(type,'Chemotherapy','Immunotherapy','Targeted therapy',
                            'Detoxificant','Experimental therapy','Hormonal therapy','Nuclear therapy'))
df_smc56 <- filter(df, group=='SMC56 SmV') %>%
  mutate(type=fct_relevel(type,'Chemotherapy','Immunotherapy','Targeted therapy',
                          'Detoxificant','Experimental therapy'))
fit <- survfit(Surv(time, status) ~ type, data = df_smc56)
pairwise_survdiff(Surv(time, status) ~ type, data = df_smc56) 
ggsurvplot(fit, data = df_smc56, pval = FALSE, ggtheme = theme_general, 
           palette = c('#C3A016FF','#246893FF','#0C1F4BFF','#C3D878FF','#58A787FF'),
           legend = 'bottom', legend.title = '', 
           xlab = 'Overall survival time (month)', xlim=c(0,150)) 
ggsave(paste(output_dir, '/figures/survival_treatment_overal_time_smc56smv.eps', sep = ''),
       dpi = 320, width = 5.5, height = 4, units = 'in')

## SmV vs others - treatment specific ----
fit <- survdiff(Surv(time, status) ~ group, data=df)
ggsurvplot_facet(fit, data = as.data.frame(df), pval = TRUE,
                 facet.by = 'type',palette = c('#A52A2A','#28446f'),
                 xlab = 'Overall survival time (months)')
ggsave(paste(output_dir, '/figures/survival_treatment_smv_others_overal_time.eps', sep = ''), 
       dpi = 320)
## immunotherapy, high vs low TMB ----
df <- inner_join(condense_survival, post_biopsy, 
                 by = join_by(sampleId), relationship = 'one-to-many') %>%
  inner_join(tmb_master, by=join_by(sampleId==Tumor_Sample_Barcode)) %>%
  select(all_of(c('patientId','time','status','type','TMB'))) %>%
  unique() %>% 
  filter(!is.na(time)) %>% 
  filter(type=='Immunotherapy') %>%
  mutate(group=with(., ifelse(TMB < 100, 'low mutational burden', 'hypermutated'))) 
fit <- survfit(Surv(time, status) ~ group, data = df)
ggsurvplot(fit, data = df, pval = T, ggtheme = theme_general, 
           legend = 'bottom', legend.title = '', 
           xlab = 'Overall survival time (month)', xlim=c(0,100)) 
ggsave(paste(output_dir, '/figures/survival_immunotherapy_low_vs_highTMB_overal_time.eps', sep = ''),
       dpi = 320, width = 5.5, height = 4, units = 'in')
## immunotherapy low vs high TMB, and SMC5/6 SmV vs others ----
df <- inner_join(condense_survival, post_biopsy, 
                 by = join_by(sampleId), relationship = 'one-to-many') %>%
  inner_join(tmb_master, by=join_by(sampleId==Tumor_Sample_Barcode)) %>%
  filter(!is.na(time)) %>% 
  filter(type=='Immunotherapy') %>%
  mutate(class=with(., ifelse(TMB >= 100,
                              ifelse(sampleId %in% smc56_smv_tumors$Tumor_Sample_Barcode,
                                     'SMC5/6 SmV + hypermutation','others + hypermutation'),
                              ifelse(sampleId %in% smc56_smv_tumors$Tumor_Sample_Barcode,
                                     'SMC5/6 SmV + low mutational burden','others + low mutational burden')))) %>%
  select(all_of(c('patientId','time','status','class'))) %>%
  unique() 
fit <- survfit(Surv(time, status) ~ class, data = df)
pairwise_survdiff(Surv(time, status) ~ class, data = df) 
ggsurvplot(fit, data = df, pval = T, ggtheme = theme_general, 
           legend = 'bottom', legend.title = '', 
           xlab = 'Overall survival time (month)', xlim=c(0,100)) 
ggsave(paste(output_dir, '/figures/survival_immunotherapy_low_vs_highTMB_overal_time_smc56smv_others.eps', sep = ''),
       dpi = 320, width = 5.5, height = 4, units = 'in')
# POST-BIOPSY survival time ----
# conduct a new survival dataset
post_biopsy_survival <- inner_join(clinical[,c('sampleId','hmfPatientId','gender','deathDate','biopsyDate',
                                          'hasSystemicPreTreatment','hasRadiotherapyPreTreatment')],
                              treatment[,c('sampleId','startDate','responseDate')],
                              by=join_by(sampleId), relationship='one-to-many') %>%
  mutate(status = with(., ifelse(is.na(deathDate),0,1))) %>%
  group_by(sampleId) %>%
  mutate(lastResponseDate = max(responseDate)) %>%
  ungroup() %>%
  mutate(survival_time = with(., ifelse(is.na(deathDate),
                                        as.duration(startDate %--% lastResponseDate)/dmonths(1),
                                        as.duration(startDate %--% deathDate)/dmonths(1)))) %>%
  group_by(sampleId) %>%
  mutate(time = max(survival_time)) %>%
  ungroup() %>%
  select(all_of(c('sampleId','time','status'))) %>%
  unique()
## treatment only ----
df <- inner_join(post_biopsy_survival, post_biopsy, 
                 by = join_by(sampleId), relationship = 'one-to-many') %>%
  select(all_of(c('patientId','time','status','type'))) %>%
  unique() %>% 
  filter(!is.na(time))
fit <- survfit(Surv(time, status) ~ type, data = df)
pairwise_survdiff(Surv(time, status) ~ type, data = df) 
ggsurvplot(fit, data = df, pval = FALSE, ggtheme = theme_general,
           palette = c('#C3A016FF','#C3D878FF','#58A787FF','#8EBACDFF','#246893FF','#163274FF','#0C1F4BFF'),
           legend = 'bottom', legend.title = '',
           xlab = 'Post biopsy survival time (month)', xlim=c(0,100)) 
ggsave(paste(output_dir, '/figures/survival_treatment_post_biopsy_time.eps', sep = ''),
       dpi = 320, width = 5.5, height = 4, units = 'in')
## treatment only, with SMC5/6 SmVs ----
df <- inner_join(post_biopsy_survival, post_biopsy, by=join_by(sampleId), 
                 relationship = 'one-to-many') %>%
  filter(!is.na(time)) %>%
  unique() %>% 
  mutate(group = with(., ifelse(sampleId %in% smc56_smv_tumors$Tumor_Sample_Barcode,'SMC56 SmV','others'))) %>%
  select(all_of(c('patientId','time','status','type','group'))) %>%
  unique() %>% 
  filter(!is.na(type)) %>%
  mutate(group=fct_relevel(group,'SMC56 SmV','others')) %>%
  mutate(type = fct_relevel(type,'Chemotherapy','Immunotherapy','Targeted therapy',
                            'Detoxificant','Experimental therapy','Hormonal therapy','Nuclear therapy'))
df_smc56 <- filter(df, group=='SMC56 SmV') %>%
  mutate(type=fct_relevel(type,'Chemotherapy','Immunotherapy','Targeted therapy',
                          'Detoxificant','Experimental therapy'))
fit <- survfit(Surv(time, status) ~ type, data = df_smc56)
pairwise_survdiff(Surv(time, status) ~ type, data = df_smc56) 
ggsurvplot(fit, data = df_smc56, pval = FALSE, ggtheme = theme_general, 
           palette = c('#C3A016FF','#246893FF','#0C1F4BFF','#C3D878FF','#58A787FF'),
           legend = 'bottom', legend.title = '', 
           xlab = 'Post biopsy survival time (month)', xlim=c(0,100)) 
ggsave(paste(output_dir, '/figures/survival_treatment_post_biopsy_time_smc56smv.eps', sep = ''),
       dpi = 320)

## SmV vs others ----
fit <- survdiff(Surv(time, status) ~ group, data=df)
pairwise_survdiff(Surv(time, status) ~ group, data = df)
ggsurvplot_facet(fit, data = as.data.frame(df), pval = TRUE,
                 facet.by = 'type',palette = c('#A52A2A','#28446f'),
                 xlab = 'Post biopsy survival time (months)')
ggsave(paste(output_dir, '/figures/survival_treatment_smv_others_post_biopsy_time.eps', sep = ''),
       dpi = 320)