# August 21st, 2025
# TMB of Hartwig tumors 

# ENVIRONMENT ####
library('tidyverse')
library('ggplot2')
library('gtsummary')
library('ggpubr',
        lib.loc = '/storage2/fs1/abby.green/Active/R_libraries/Bunny-Wunnies Freak Out.Bunny-Wunnies Freak Out')
library('VennDiagram')
theme_general <- theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

input_dir <- '/storage2/fs1/abby.green/Active/Hartwig'
output_dir <- '/storage2/fs1/abby.green/Active/Users/thi/SMC56_de_novo/2025/Hartwig'
smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
                'NSMCE4A', 'NSMCE4B', 'EID3',
                'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')

# Calculate TMB and output it ----
smv_master <- list.files(paste(output_dir,'/data',sep=''), 
                         full.names=T, pattern='patients_pass_only') %>%
  map(read.csv, header=T) %>%
  bind_rows()
# based on the smv_master file
tmb_master <- filter(smv_master, FILTER=='PASS') %>%
  select(all_of(c('Tumor_Sample_Barcode','Chromosome','Start_Position',
                  'End_Position','Reference_Allele','Tumor_Seq_Allele2'))) %>%
  unique() %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(abs_mut_number=n()) %>%
  mutate(TMB=abs_mut_number/2800)
write.csv(tmb_master, paste(output_dir, '/data/pass_tmb_master.csv', sep = ''),
          row.names = F)
# INPUT ####
smv_master <- list.files(paste(output_dir,'/data',sep=''), 
                         full.names=T, pattern='patients_pass_only') %>%
  map(read.csv, header=T) %>%
  bind_rows()
tmb_master <- read.csv(paste(output_dir,'/data/pass_tmb_master.csv',sep=''),header=T)
post_biopsy <- read_tsv(paste(input_dir, '/post_biopsy_drugs.tsv', sep = ''), na = c('null', ''))
# load the cohort files generated on Oct 28th
cohort_master <- list.files(paste(output_dir, '/data', sep = ''), full.names=T, pattern = 'cohort.csv$') %>%
  map(read.csv, header = T) %>%
  bind_rows()
# trinucleotide master 
trinu_master <- read.csv(paste(output_dir, '/data/indv_trinu_context.csv', sep = ''),
                               header = TRUE)
# just to make sure that I dont take into account 
# patients without survival information 
condense_survival <- read.csv(paste(output_dir, '/data/survival.csv', sep = ''),
                              header = TRUE) %>% 
  filter(!is.na(time))
# SMC5/6 alterations ----
tmb_cohort <- inner_join(tmb_master, cohort_master, 
                      by=join_by(Tumor_Sample_Barcode), relationship = 'one-to-many') %>%
  filter(Tumor_Sample_Barcode %in% condense_survival$sampleId) %>% 
  select(all_of(c('Tumor_Sample_Barcode','TMB','COHORT'))) %>%
  unique() %>%
  mutate(color=with(.,ifelse(COHORT %in% c('del','nsnd','syn'),
                             'SmV',ifelse(COHORT %in% c('high_cnv','low_cnv'),
                                          'CNV','non-altered')))) %>%
  select(all_of(c('Tumor_Sample_Barcode','color','TMB'))) %>%
  unique() %>%
  mutate(color=fct_relevel(color,'SmV','CNV','non-altered')) 
ggplot(tmb_cohort,aes(x=color, y=log10(TMB))) +
  geom_boxplot(aes(color=color)) +
  theme_general +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#A52A2A','#702963','#28446f')) +
  stat_compare_means(method = 'wilcox', label = 'p.format',
                     comparisons = list(c('SmV','CNV'),c('CNV','non-altered'),c('SmV','non-altered')))
ggsave(paste(output_dir, '/figures/log10tmb_3groups.eps', sep = ''),
       dpi = 320, width = 4, height = 4, units = 'in')
# SMC5/6 intronic vs exonic SmV ----
smc56_exon_smv <- filter(smv_master,
                         Hugo_Symbol %in% smc56_gene & Variant_Classification != 'Intron') %>%
  mutate(COHORT='SMC5/6 exon mutations')
smc56_intron_smv <- filter(smv_master, 
                           Hugo_Symbol %in% smc56_gene & Variant_Classification == 'Intron') %>%
  filter(!Tumor_Sample_Barcode %in% smc56_exon_smv$Tumor_Sample_Barcode) %>%
  mutate(COHORT='SMC5/6 intronic mutations')
tmb_intron_exon <- bind_rows(unique(smc56_exon_smv[,c('Tumor_Sample_Barcode','COHORT')]),
                             unique(smc56_intron_smv[,c('Tumor_Sample_Barcode','COHORT')]),
                             unique(cohort_master[,c('Tumor_Sample_Barcode','COHORT')])) %>%
  filter(COHORT %in% c('SMC5/6 intronic mutations','SMC5/6 exon mutations',
                       'low_cnv','high_cnv','wt')) %>%
  inner_join(tmb_master, by=join_by(Tumor_Sample_Barcode), relationship='many-to-one') %>%
  filter(Tumor_Sample_Barcode %in% condense_survival$sampleId) %>%
  mutate(color=with(.,ifelse(COHORT=='wt','non-altered',
                             ifelse(COHORT %in% c('low_cnv','high_cnv'),'CNV',COHORT)))) %>%
  select(all_of(c('Tumor_Sample_Barcode','color','TMB'))) %>%
  unique() %>%
  mutate(color=fct_relevel(color,'SMC5/6 intronic mutations','SMC5/6 exon mutations',
                           'CNV','non-altered'))
ggplot(tmb_intron_exon,aes(x=color,y=log10(TMB))) +
  geom_boxplot() +
  theme_general +
  theme(legend.position = 'none') +
  stat_compare_means(comparisons=list(c('SMC5/6 intronic mutations','SMC5/6 exon mutations'),
                                      c('SMC5/6 exon mutations','CNV'),c('CNV','non-altered'),
                                      c('SMC5/6 intronic mutations','CNV'),
                                      c('SMC5/6 exon mutations','non-altered'),
                                      c('SMC5/6 intronic mutations','non-altered')),
                     method = 'wilcox', label = 'p.format')
ggsave(paste(output_dir,'/figures/log10TMB_intron_exonSMC56.eps',sep=''), 
       dpi=320, width=4, height=4, units='in')
# different post-biopsy treatment -----
tmb_treatment <- inner_join(tmb_master, post_biopsy,
                      by=join_by(Tumor_Sample_Barcode==sampleId), relationship='one-to-many') %>%
  filter(!is.na(type) & Tumor_Sample_Barcode %in% condense_survival$sampleId) %>%
  select(all_of(c('Tumor_Sample_Barcode','TMB','type'))) %>%
  unique() 
ggplot(tmb_treatment,aes(x=type, y=log10(TMB), color=type)) +
  geom_boxplot() +
  theme_general +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#C3A016FF','#C3D878FF','#58A787FF','#8EBACDFF','#246893FF','#163274FF','#0C1F4BFF')) +
  stat_compare_means(methods = 'wilcox', comparisons = list(c('Hormonal therapy','Immunotherapy'),
                                                            c('Experimental therapy','Immunotherapy'),
                                                            c('Detoxificant','Immunotherapy'),
                                                            c('Chemotherapy','Immunotherapy'),
                                                            c('Nuclear therapy','Immunotherapy'),
                                                            c('Targeted therapy','Immunotherapy')),
                     hide.ns = TRUE) 
ggsave(paste(output_dir, '/figures/log10tmb_treatment.eps', sep = ''),
       dpi = 320, width = 7, height = 4, units = 'in')
# treatment but into SmV and others ----
tmb_cohort_treatment <- inner_join(tmb_master, post_biopsy, 
                      by = join_by(Tumor_Sample_Barcode==sampleId), relationship = 'one-to-many') %>%
  inner_join(cohort_master, by=join_by(Tumor_Sample_Barcode), relationship='many-to-many') %>%
  filter(!is.na(type) & Tumor_Sample_Barcode %in% condense_survival$sampleId) %>%
  mutate(color=with(.,ifelse(COHORT %in% c('del','nsnd','syn'),'SmV','others'))) %>%
  select(all_of(c('Tumor_Sample_Barcode','color','type','TMB'))) %>%
  unique() %>%
  mutate(color=fct_relevel(color,'SmV','others')) 
ggplot(tmb_cohort_treatment, aes(x=type, y=log10(TMB), color=color)) +
  geom_boxplot() +
  theme_general +
  scale_color_manual(values = c('#A52A2A','#28446f')) +
  stat_compare_means(method = 'wilcox', hide.ns = TRUE, label='p.format')
ggsave(paste(output_dir, '/figures/log10tmb_treatment_smv_others.eps', sep = ''),
       dpi = 320, width = 7, height = 4, units = 'in')
