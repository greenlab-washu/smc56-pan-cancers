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
library('ggpubr')
library('survminer')
theme_general <- theme_classic() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())
# this directory contains secondary data derived from raw data provided by Hartwig
# whose generation was documented in data_cleaning.R
# the secondary data does not come with this GitHub page 
input_dir <- 'Hartwig' 
output_dir <- 'Hartwig'

smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
                'NSMCE4A', 'NSMCE4B', 'EID3',
                'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')

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
# based on post-biopsy treatment
# OVERAL survival time ----

## Fig4d - treatment only, with SMC5/6 SmVs ----
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
## Fig 4e - SmV vs others - treatment specific ----
fit <- survdiff(Surv(time, status) ~ group, data=df)
ggsurvplot_facet(fit, data = as.data.frame(df), pval = TRUE,
                 facet.by = 'type',palette = c('#A52A2A','#28446f'),
                 xlab = 'Overall survival time (months)')
