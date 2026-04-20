# ENVIRONMENT ----
library('tidyverse')
library('ggplot2')
library('ggVennDiagram')
library('ggpubr')
library('survminer')
library('survival')
library('MutationalPatterns')
theme_general <- theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = 'bottom')

input_dir <- 'TCGA/data'
output_dir <- 'TCGA'
smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
                'NSMCE4A', 'NSMCE4B', 'EID3',
                'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')
mmr_genes <- c('MLH1','MLH2','MLH3','MSH2','MSH6','PMS1','PMS2')
pole_genes <- c('POLE','POLE2','POLE3','POLE4')
cosmic_v3.2 <- get_known_signatures(muttype='snv', source='COSMIC_v3.2', genome='GRCh38')

color_code <- c('#FFC72CFF','#FF9E1BFF',
                '#E87722FF','#DC582AFF',
                '#EEF3FFFF','#CDE5F9FF',
                '#85B6CEFF',
                '#64A8A8FF','#4A9152FF','#45681EFF','#54450FFF','#362904FF',
                '#707271FF',
                c(rep('#dddddc',47)))
names(color_code) <- c('SBS2','SBS13', #apobec
                       'SBS84','SBS85', #AID
                       'SBS10a', 'SBS10b', #PolE mutation
                       'SBS14', #PolE+MMRd
                       'SBS6','SBS15','SBS21','SBS26','SBS44', #MMRd
                       'SBS87', #thiopurine chemotherapy 
                       'SBS1','SBS5', #Clock-like
                       'SBS4','SBS29','SBS92', #tobaco smoking
                       'SBS11','SBS25','SBS31','SBS35','SBS86',#chemotherapy
                       "SBS3","SBS7a","SBS7b","SBS7c","SBS7d","SBS8",
                       "SBS9","SBS10c","SBS10d","SBS12","SBS16","SBS17a",
                       "SBS17b","SBS18","SBS19","SBS20","SBS22","SBS23",
                       "SBS24","SBS28","SBS30","SBS32","SBS33","SBS34",
                       "SBS36","SBS37","SBS38","SBS39","SBS40","SBS41",
                       "SBS42","SBS88","SBS89","SBS90","SBS91","SBS93","SBS94")

# INPUT ----
snv_master <- read.csv(paste(input_dir, '/gdc_snv.csv', sep = ''), header = TRUE)
tmb_master <- read.csv(paste(output_dir, '/data/gdc_smv30.csv', sep = ''),header=T)
trinu_master <- read.csv(paste(input_dir, '/master_96.csv', sep = ''), header=T)
colnames(trinu_master) <- str_replace_all(colnames(trinu_master), '\\.', '-')
# SMC5/6 cohorts
cohort_sample <- 
  list.files(input_dir, full.names = TRUE, pattern = 'cohort_for_fpkm.csv$') %>%
  purrr::map(read.csv, header = TRUE) 
colnames(cohort_sample[[2]]) <- colnames(cohort_sample[[4]])
cohort_master <- bind_rows(cohort_sample)
# patients clinical files - survival
tcga_survival <- list.files(path = dir('cBioPortal/pan_can_atlas',
                                       full.names = TRUE),
                            full.names = TRUE, pattern = 'data_clinical_patient.txt') %>%
  map(read.delim, sep = '\t', header = TRUE, comment.char = '#', fill = FALSE) %>%
  purrr::map(mutate_at, c('DAYS_LAST_FOLLOWUP', 'DAYS_TO_BIRTH', 'AJCC_STAGING_EDITION'), as.integer) %>%
  bind_rows()

# 5b - Venn diagram ----
smc56_tumors <- filter(cohort_master, 
                       Hugo_Symbol %in% smc56_gene & alteration_type %in% c('nonsyn_nondel_mutation','del_mutation')) %>%
  mutate(sample_id=str_extract(Tumor_Sample_Barcode,'^.{15}'))
mmr_tumors <- filter(cohort_master,
                     Hugo_Symbol %in% mmr_genes & alteration_type %in% c('nonsyn_nondel_mutation','del_mutation')) %>%
  mutate(sample_id=str_extract(Tumor_Sample_Barcode,'^.{15}'))
pole_tumors <- filter(cohort_master,
                      Hugo_Symbol %in% pole_genes & alteration_type %in% c('nonsyn_nondel_mutation','del_mutation')) %>%
  mutate(sample_id=str_extract(Tumor_Sample_Barcode,'^.{15}'))
x <- list(
  smc56_nonsyn=unique(smc56_tumors$sample_id),
  mmr_nonsyn=unique(mmr_tumors$sample_id),
  pole_nonsyn=unique(pole_tumors$sample_id)
)
ggVennDiagram(x, label_alpha = 0, label_font = 'Helvetica') + # need to specify this font, or else I could not save the VennDiagram wtih ggsave
  scale_fill_gradient(low='white',high='#A52A2A')
ggsave(paste(output_dir,'/figures/venn_pole_mmr_smc56_nonsyn_v2.eps', sep=''),
       dpi = 320, width = 5, height = 5, units = 'in')
# 5c - log10TMB - SMC5/6 nonsyn tumors ----
df <- filter(tmb_master, str_detect(Tumor_Sample_Barcode,'A$')==T) %>%
  mutate(sample_id=str_extract(Tumor_Sample_Barcode,'^.{15}')) %>%
  filter(sample_id %in% smc56_tumors$sample_id) %>%
  mutate(class=with(.,ifelse(sample_id %in% intersect(pole_tumors$sample_id, mmr_tumors$sample_id),
                             'SMC56, POLE, and MMR nonsyn mutations',
                             ifelse(sample_id %in% setdiff(mmr_tumors$sample_id, pole_tumors$sample_id),
                                    'SMC56 and MMR nonsyn mutations',
                                    ifelse(sample_id %in% setdiff(pole_tumors$sample_id, mmr_tumors$sample_id),
                                           'SMC56 and POLE nonsyn mutations',
                                           'SMC56 nonsyn mutations only'))))) %>%
  select(all_of(c('sample_id','TMB','class'))) %>%
  unique() %>% 
  mutate(class=fct_relevel(class,'SMC56 nonsyn mutations only',
                           'SMC56 and MMR nonsyn mutations','SMC56 and POLE nonsyn mutations',
                           'SMC56, POLE, and MMR nonsyn mutations'))
ggplot(df, aes(x=class, y=log10(TMB), color=class)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c('SMC56 nonsyn mutations only','SMC56 and MMR nonsyn mutations'),
                                        c('SMC56 nonsyn mutations only','SMC56 and POLE nonsyn mutations'),
                                        c('SMC56 nonsyn mutations only','SMC56, POLE, and MMR nonsyn mutations'),
                                        c('SMC56 and MMR nonsyn mutations','SMC56 and POLE nonsyn mutations'),
                                        c('SMC56 and MMR nonsyn mutations','SMC56, POLE, and MMR nonsyn mutations'),
                                        c('SMC56 and POLE nonsyn mutations','SMC56, POLE, and MMR nonsyn mutations')),
                     method='wilcox', label='p.adj') +
  theme_general +
  scale_color_manual(values = c('#FDE725FF','#7AD151FF','#22A884FF','#2A788EFF'))
ggsave(paste(output_dir, '/figures/log10TMB_pole_mmr_smc56_nonsyn_v2.eps', sep = ''),
       dpi = 320, width = 5, height = 5, units = 'in')  
# 5d - Survival ----
df <- mutate(df, patient_id=str_extract(sample_id, '^.{12}'))
df_survival <- filter(tcga_survival, !is.na(OS_MONTHS)) %>%
  inner_join(df, by=join_by(PATIENT_ID==patient_id)) %>%
  separate(OS_STATUS, c('OS_STATUS_DIGIT','OS_STATUS'), sep=':') %>%
  mutate_at(c('OS_STATUS_DIGIT'), as.numeric) %>%
  select(all_of(c('PATIENT_ID','class','OS_MONTHS','OS_STATUS_DIGIT'))) %>%
  unique()
fit <- df_survival %>%
  surv_fit(Surv(OS_MONTHS, OS_STATUS_DIGIT) ~ class, data = .)
df_survival %>% pairwise_survdiff(Surv(OS_MONTHS, OS_STATUS_DIGIT) ~ class, data = .)
df_survival %>%
  ggsurvplot(fit, data = ., pval=T, # there is no statical significance here so don't need to show p value
             palette = c('#FDE725FF','#7AD151FF','#22A884FF','#2A788EFF'), 
             ggtheme = theme_general,
             legend = 'right', legend.title = '',
             xlab = 'Time (month)') 
ggsave(paste(output_dir, '/figures/survival_pole_mmr_smc56_nonsyn_v2.eps', sep=''),
       dpi = 320, width = 7, height = 5, units = 'in')
# S9b - pool SBSsig ---- 
df1 <- select(trinu_master, all_of(c(setdiff(smc56_tumors$Tumor_Sample_Barcode, # this is SMC5/6 only 
                                             union(mmr_tumors$Tumor_Sample_Barcode, 
                                                   pole_tumors$Tumor_Sample_Barcode)))))
df2 <- select(trinu_master, all_of(c(intersect(smc56_tumors$Tumor_Sample_Barcode, # this is SMC5/6 + MMR genes
                                               setdiff(mmr_tumors$Tumor_Sample_Barcode, 
                                                       pole_tumors$Tumor_Sample_Barcode)))))
df3 <- select(trinu_master, all_of(c(intersect(smc56_tumors$Tumor_Sample_Barcode, # this is SMC5/6 + POLE 
                                               setdiff(pole_tumors$Tumor_Sample_Barcode, 
                                                       mmr_tumors$Tumor_Sample_Barcode)))))
df4 <- select(trinu_master, all_of(c(intersect(smc56_tumors$Tumor_Sample_Barcode, # this is SMC5/6 + MMR + POLE
                                               intersect(mmr_tumors$Tumor_Sample_Barcode, 
                                                         pole_tumors$Tumor_Sample_Barcode))))) 
pool_trinu <- map(list(df1, df2, df3, df4), rowSums) %>%
  map(as.data.frame) %>%
  bind_cols()
colnames(pool_trinu) <- c('df1','df2','df3','df4')
fit <- fit_to_signatures(pool_trinu, cosmic_v3.2)
df_plot <- apply(fit$contribution, 2, function(x) x / sum(x)) %>%
  as_tibble(rownames = 'SBSsig') %>%
  pivot_longer(!c('SBSsig'), names_to = 'class', values_to = 'rel_cont') %>%
  mutate(SBSsig = fct_relevel(SBSsig, 'SBS2','SBS13', #apobec
                              'SBS84','SBS85', #AID
                              'SBS10a', 'SBS10b', #PolE mutation
                              'SBS14', #PolE+MMRd
                              'SBS6','SBS15','SBS21','SBS26','SBS44', #MMRd
                              'SBS87', #thiopurine chemotherapy 
                              'SBS1','SBS5', #Clock-like
                              'SBS4','SBS29','SBS92', #tobaco smoking
                              'SBS11','SBS25','SBS31','SBS35','SBS86',#chemotherapy
                              "SBS3","SBS7a","SBS7b","SBS7c","SBS7d","SBS8",
                              "SBS9","SBS10c","SBS10d","SBS12","SBS16","SBS17a",
                              "SBS17b","SBS18","SBS19","SBS20","SBS22","SBS23",
                              "SBS24","SBS28","SBS30","SBS32","SBS33","SBS34",
                              "SBS36","SBS37","SBS38","SBS39","SBS40","SBS41",
                              "SBS42","SBS88","SBS89","SBS90","SBS91","SBS93","SBS94")
  )
ggplot(df_plot, aes(x=class, y=rel_cont, fill=SBSsig)) +
  geom_col(position = 'stack') +
  theme_general +
  scale_fill_manual(values = color_code) +
  theme(legend.position = 'right')
ggsave(paste(output_dir, '/figures/pool_sbssig_pole_mmr_smc56_nonsyn_v2.eps', sep = ''),
       dpi = 320, width = 7, height = 7, units = 'in')
# 5e - survival - POLE + MMR genes nonsyn mutations in SMC5/6 non-altered tumors -----
wt_smc56_tumors <- filter(cohort_master,COHORT=='wt') %>%
  mutate(patient_id=str_extract(Tumor_Sample_Barcode, '^.{12}'))
mmr_tumors <- mutate(mmr_tumors, patient_id=str_extract(Tumor_Sample_Barcode, '^.{12}'))
pole_tumors <- mutate(pole_tumors, patient_id=str_extract(Tumor_Sample_Barcode, '^.{12}'))
df_survival <- filter(tcga_survival, !is.na(OS_MONTHS) & PATIENT_ID %in% wt_smc56_tumors$patient_id) %>%
  separate(OS_STATUS, c('OS_STATUS_DIGIT','OS_STATUS'), sep=':') %>%
  mutate_at(c('OS_STATUS_DIGIT'), as.numeric) %>%
  mutate(class=with(., ifelse(PATIENT_ID %in% intersect(mmr_tumors$patient_id, pole_tumors$patient_id),
                              'MMR and POLE nonsyn mutations', ifelse(PATIENT_ID %in% setdiff(mmr_tumors$patient_id, pole_tumors$patient_id),
                                                                      'MMR nonsyn mutations only', ifelse(PATIENT_ID %in% setdiff(pole_tumors$patient_id, mmr_tumors$patient_id),
                                                                                                          'POLE nonsyn mutations only', 'non mutated'))))) %>%
  select(all_of(c('PATIENT_ID','class','OS_MONTHS','OS_STATUS_DIGIT'))) %>%
  unique() %>% 
  mutate(class=fct_relevel(class, 'non mutated','MMR nonsyn mutations only',
                           'POLE nonsyn mutations only','MMR and POLE nonsyn mutations'))
fit <- df_survival %>%
  surv_fit(Surv(OS_MONTHS, OS_STATUS_DIGIT) ~ class, data = .)
df_survival %>% pairwise_survdiff(Surv(OS_MONTHS, OS_STATUS_DIGIT) ~ class, data = .)
df_survival %>%
  ggsurvplot(fit, data = ., pval=T, # there is no statical significance here so don't need to show p value
             palette = c('#AFD0C9FF','#91A1BAFF','#51628EFF','#182044FF'), 
             ggtheme = theme_general,
             legend = 'right', legend.title = '',
             xlab = 'Time (month)', xlim = c(0,250)) 
ggsave(paste(output_dir, '/figures/survival_pole_mmr_nonsyn_wtSMC56_v2.eps', sep=''),
       dpi = 320, width = 7, height = 5, units = 'in')
# S9d - survival - POLE mutations + SMC5/6 non-altered vs SmV ----
smc56_tumors <- filter(cohort_master, COHORT %in% c('syn','del','nonsyn_nondel')) %>%
  mutate(patient_id=str_extract(Tumor_Sample_Barcode,'^.{12}'))
df_survival <- filter(tcga_survival, !is.na(OS_MONTHS) & PATIENT_ID %in% union(smc56_tumors$patient_id, wt_smc56_tumors$patient_id)) %>%
  separate(OS_STATUS, c('OS_STATUS_DIGIT','OS_STATUS'), sep=':') %>%
  mutate_at(c('OS_STATUS_DIGIT'), as.numeric) %>%
  mutate(class=with(., ifelse(PATIENT_ID %in% pole_tumors$patient_id,
                              ifelse(PATIENT_ID %in% smc56_tumors$patient_id,'SMC5/6 SmV + POLE nonsyn mutations','SMC5/6 non-altered + POLE nonsyn mutations'),
                              ifelse(PATIENT_ID %in% smc56_tumors$patient_id,'SMC5/6 SmV + no POLE nonsyn mutations','SMC5/6 non-altered + no POLE nonsyn mutations')))) %>%
  select(all_of(c('PATIENT_ID','class','OS_MONTHS','OS_STATUS_DIGIT'))) %>%
  unique() %>% 
  mutate(class=fct_relevel(class, 'SMC5/6 non-altered + POLE nonsyn mutations',
                           'SMC5/6 SmV + POLE nonsyn mutations',
                           'SMC5/6 non-altered + no POLE nonsyn mutations',
                           'SMC5/6 SmV + no POLE nonsyn mutations'))
fit <- df_survival %>%
  surv_fit(Surv(OS_MONTHS, OS_STATUS_DIGIT) ~ class, data = .)
df_survival %>% pairwise_survdiff(Surv(OS_MONTHS, OS_STATUS_DIGIT) ~ class, data = .)
df_survival %>%
  ggsurvplot(fit, data = ., pval=T, # there is no statical significance here so don't need to show p value
             palette = c('#088BBEFF','#172869FF','#F28A8AFF','#D9565CFF'), 
             ggtheme = theme_general,
             legend = 'right', legend.title = '',
             xlab = 'Overall survival time (month)', xlim = c(0,250)) 
ggsave(paste(output_dir, '/figures/survival_pole_smc56smv_wt_v2.eps', sep=''),
       dpi = 320, width = 7, height = 5, units = 'in')
# S9c - survival - MMR genes mutations + SMC5/6 non-altered vs SmV ----
df_survival <- filter(tcga_survival, !is.na(OS_MONTHS) & PATIENT_ID %in% union(smc56_tumors$patient_id, wt_smc56_tumors$patient_id)) %>%
  separate(OS_STATUS, c('OS_STATUS_DIGIT','OS_STATUS'), sep=':') %>%
  mutate_at(c('OS_STATUS_DIGIT'), as.numeric) %>%
  mutate(class=with(., ifelse(PATIENT_ID %in% mmr_tumors$patient_id,
                              ifelse(PATIENT_ID %in% smc56_tumors$patient_id,'SMC5/6 SmV + MMR nonsyn mutations','SMC5/6 non-altered + MMR nonsyn mutations'),
                              ifelse(PATIENT_ID %in% smc56_tumors$patient_id,'SMC5/6 SmV + no MMR nonsyn mutations','SMC5/6 non-altered + no MMR nonsyn mutations')))) %>%
  select(all_of(c('PATIENT_ID','class','OS_MONTHS','OS_STATUS_DIGIT'))) %>%
  unique() %>% 
  mutate(class=fct_relevel(class, 'SMC5/6 non-altered + MMR nonsyn mutations',
                           'SMC5/6 SmV + MMR nonsyn mutations',
                           'SMC5/6 non-altered + no MMR nonsyn mutations',
                           'SMC5/6 SmV + no MMR nonsyn mutations'))
fit <- df_survival %>%
  surv_fit(Surv(OS_MONTHS, OS_STATUS_DIGIT) ~ class, data = .)
df_survival %>% pairwise_survdiff(Surv(OS_MONTHS, OS_STATUS_DIGIT) ~ class, data = .)
df_survival %>%
  ggsurvplot(fit, data = ., pval=T, # there is no statical significance here so don't need to show p value
             palette = c('#088BBEFF','#172869FF','#F28A8AFF','#D9565CFF'), 
             ggtheme = theme_general,
             legend = 'right', legend.title = '',
             xlab = 'Overall survival time (month)', xlim = c(0,250)) 
ggsave(paste(output_dir, '/figures/survival_mmr_smc56smv_wt.eps', sep=''),
       dpi = 320, width = 7, height = 5, units = 'in')
