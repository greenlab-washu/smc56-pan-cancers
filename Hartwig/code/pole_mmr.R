# December 12th, 2025 
# SBSsig of patients with nonsyn MMR and POLE genes mutations
# all patients and only those receiving immunotherapy

# ENVIRONMENT ####
library('MutationalPatterns')
library('BSgenome')
library('tidyverse')
library('ggplot2')
library('lubridate')
library('survival')
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
#load reference genome - GRCh37, which is the reference genome used in Hartwig
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg19' 
library(ref_genome, character.only = TRUE)
mmr_genes <- c('MLH1','MLH2','MLH3','MSH2','MSH6','PMS1','PMS2')
pole_genes <- c('POLE','POLE2','POLE3','POLE4')
# load COSMIC signature
snv_cosmic_3.2 <- get_known_signatures(muttype = 'snv', 'COSMIC_v3.2', genome = 'GRCh37')
indel_cosmic_3.2 <- get_known_signatures(muttype = 'indel', 'COSMIC_v3.2', genome = 'GRCh37')
## color ----
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
clinical <- read_tsv(paste(input_dir, '/metadata.tsv', sep = ''), na = c('null', ''))
post_biopsy <- read_tsv(paste(input_dir, '/post_biopsy_drugs.tsv', sep = ''), na = c('null', ''))
treatment <- read_tsv(paste(input_dir, '/treatment_responses.tsv', sep = ''), na = c('null', '')) 
# just to make sure that I dont take into account 
# patients without survival information 
condense_survival <- read.csv(paste(output_dir, '/data/survival.csv', sep = ''),
                              header = TRUE) %>% 
  filter(!is.na(time))
# load mutational matrix of indv patients 
mutation_matrix <- read.csv(paste(output_dir,'/data/indv_trinu_context.csv',sep=''),header = T)
# all PASS SmV in Hartwig
pass_smv_master <- list.files(paste(output_dir, '/data', sep = ''), full.names = T, pattern = 'patients_pass_only') %>%
  map(read.csv, header = T) %>%
  bind_rows()
# SBSsig ----
mmr_tumors <- filter(pass_smv_master, Hugo_Symbol %in% mmr_genes & 
                       Tumor_Sample_Barcode %in% condense_survival$sampleId & # only samples that passed QC
                       HGVSp_Short != '' & str_detect(HGVSp_Short, '=')==F) # missense mutations only
pole_tumors <- filter(pass_smv_master, Hugo_Symbol %in% pole_genes &
                        Tumor_Sample_Barcode %in% condense_survival$sampleId & # only samples that passed QC
                        HGVSp_Short != '' & str_detect(HGVSp_Short,'=')==F) # missense mutations only 
df1 <- select(mutation_matrix, all_of(setdiff(condense_survival$sampleId,
                                              union(mmr_tumors$Tumor_Sample_Barcode,
                                                    pole_tumors$Tumor_Sample_Barcode)))) # neither MMR nor POLE genes mutations
df2 <- select(mutation_matrix, all_of(setdiff(mmr_tumors$Tumor_Sample_Barcode,
                                              pole_tumors$Tumor_Sample_Barcode))) # MMR genes mutations only
df3 <- select(mutation_matrix, all_of(setdiff(pole_tumors$Tumor_Sample_Barcode,
                                              mmr_tumors$Tumor_Sample_Barcode))) # POLE genes mutations only
df4 <- select(mutation_matrix, all_of(intersect(mmr_tumors$Tumor_Sample_Barcode,
                                                pole_tumors$Tumor_Sample_Barcode))) # having both POLE and MMR genes mutations
pool_trinu <- map(list(df1, df2, df3, df4), rowSums) %>%
  map(as.data.frame) %>%
  bind_cols()
colnames(pool_trinu) <- c('df1','df2','df3','df4')
fit <- fit_to_signatures(pool_trinu, snv_cosmic_3.2)
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
ggsave(paste(output_dir, '/figures/pole_mmr_sbssig.eps', sep = ''), dpi=320)
# Survival ----
df <- condense_survival %>%
  mutate(pole_mmr=with(.,ifelse(sampleId %in% intersect(mmr_tumors$Tumor_Sample_Barcode,pole_tumors$Tumor_Sample_Barcode),
                                'MMR and POLE genes mutations',
                                ifelse(sampleId %in% setdiff(mmr_tumors$Tumor_Sample_Barcode,pole_tumors$Tumor_Sample_Barcode),
                                       'MMR genes mutations only',
                                       ifelse(sampleId %in% setdiff(pole_tumors$Tumor_Sample_Barcode,mmr_tumors$Tumor_Sample_Barcode),
                                              'POLE genes mutations only',
                                              'no MMR nor POLE mutations'))))) %>%
  select(all_of(c('sampleId','time','status','pole_mmr'))) %>%
  unique() %>%
  filter(!is.na(time)) %>%
  mutate(pole_mmr=fct_relevel(pole_mmr,'no MMR nor POLE mutations','MMR genes mutations only',
                              'POLE genes mutations only','MMR and POLE genes mutations'))
fit <- survfit(Surv(time, status) ~ pole_mmr, data = df)
pairwise_survdiff(Surv(time, status) ~ pole_mmr, data = df) 
ggsurvplot(fit, data = df, pval = FALSE, ggtheme = theme_general, 
           palette = c('#FDE725FF','#7AD151FF','#22A884FF','#2A788EFF'),
           legend = 'bottom', legend.title = '', 
           xlab = 'Overall survival time (month)', xlim=c(0,150)) 
ggsave(paste(output_dir, '/figures/pole_mmr_overal_survival.eps', sep = ''),
       dpi = 320, width = 5.5, height = 4, units = 'in')
# Immunotherapy patients only ----
## SBSsig ----
immuno_patients <- filter(post_biopsy, type=='Immunotherapy')
df1 <- select(mutation_matrix, all_of(setdiff(immuno_patients$sampleId,
                                              union(mmr_tumors$Tumor_Sample_Barcode,
                                                    pole_tumors$Tumor_Sample_Barcode)))) # neither MMR nor POLE genes mutations
df2 <- select(mutation_matrix, all_of(setdiff(mmr_tumors$Tumor_Sample_Barcode,
                                              pole_tumors$Tumor_Sample_Barcode))) %>% # MMR genes mutations only
  select(any_of(immuno_patients$sampleId))
df3 <- select(mutation_matrix, all_of(setdiff(pole_tumors$Tumor_Sample_Barcode,
                                              mmr_tumors$Tumor_Sample_Barcode))) %>% # POLE genes mutations only
  select(any_of(immuno_patients$sampleId))
df4 <- select(mutation_matrix, all_of(intersect(mmr_tumors$Tumor_Sample_Barcode,
                                                pole_tumors$Tumor_Sample_Barcode))) %>% # having both POLE and MMR genes mutations
  select(any_of(immuno_patients$sampleId))
pool_trinu <- map(list(df1, df2, df3, df4), rowSums) %>%
  map(as.data.frame) %>%
  bind_cols()
colnames(pool_trinu) <- c('df1','df2','df3','df4')
fit <- fit_to_signatures(pool_trinu, snv_cosmic_3.2)
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
ggsave(paste(output_dir, '/figures/immuno_pole_mmr_sbssig.eps', sep = ''), dpi=320)
## Survival ----
df <- condense_survival %>%
  mutate(pole_mmr=with(.,ifelse(sampleId %in% intersect(mmr_tumors$Tumor_Sample_Barcode,pole_tumors$Tumor_Sample_Barcode),
                                'MMR and POLE genes mutations',
                                ifelse(sampleId %in% setdiff(mmr_tumors$Tumor_Sample_Barcode,pole_tumors$Tumor_Sample_Barcode),
                                       'MMR genes mutations only',
                                       ifelse(sampleId %in% setdiff(pole_tumors$Tumor_Sample_Barcode,mmr_tumors$Tumor_Sample_Barcode),
                                              'POLE genes mutations only',
                                              'no MMR nor POLE mutations'))))) %>%
  filter(sampleId %in% immuno_patients$sampleId) %>%
  select(all_of(c('sampleId','time','status','pole_mmr'))) %>%
  unique() %>%
  filter(!is.na(time)) %>%
  mutate(pole_mmr=fct_relevel(pole_mmr,'no MMR nor POLE mutations','MMR genes mutations only',
                              'POLE genes mutations only','MMR and POLE genes mutations'))
fit <- survfit(Surv(time, status) ~ pole_mmr, data = df)
pairwise_survdiff(Surv(time, status) ~ pole_mmr, data = df) 
ggsurvplot(fit, data = df, pval = FALSE, ggtheme = theme_general, 
           palette = c('#FDE725FF','#7AD151FF','#22A884FF','#2A788EFF'),
           legend = 'bottom', legend.title = '', 
           xlab = 'Overall survival time (month)', xlim=c(0,90)) 
ggsave(paste(output_dir, '/figures/immuno_pole_mmr_overal_survival.eps', sep = ''),
       dpi = 320, width = 5.5, height = 4, units = 'in')
