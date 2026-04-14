# Sep 5th, 2025
# Just confirmed that all SmV from TCGA
# that I have access to were FILTER == PASS
# so maybe I should do the same thing for this dataset as well 


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
# INPUT ####
pass_smv_master <- list.files(paste(output_dir, '/data', sep = ''), full.names = T, pattern = 'patients_pass_only') %>%
  map(read.csv, header = T) %>%
  bind_rows()
# make sure that we only include the patients with complete clinical data 
condense_survival <- read.csv(paste(output_dir, '/data/survival.csv', sep = ''),
                              header = TRUE) %>% 
  filter(!is.na(time))
# TMB calculation -----
pass_tmb_master <- select(pass_snv_master, all_of(c(
  'Hugo_Symbol','Chromosome','Start_Position','End_Position',
  'Tumor_Sample_Barcode','Reference_Allele','Tumor_Seq_Allele2'
))) %>%
  unique() %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(abs_mut_number = n()) %>%
  mutate(TMB = abs_mut_number/2800)
write.csv(pass_tmb_master, paste(output_dir, '/data/pass_tmb_master.csv',
                                 sep = ''), row.names = FALSE)
smc56_snv <- filter(pass_snv_master, Hugo_Symbol %in% smc56_gene) %>%
  filter(Variant_Classification != 'Intron')
df_plot<- pass_tmb_master %>% 
  filter(Tumor_Sample_Barcode %in% condense_survival$sampleId) %>% # show only those with complete clinical info
  mutate(Tumor_Sample_Barcode = fct_reorder(Tumor_Sample_Barcode, -desc(TMB))) %>%
  mutate(quartile = with(., ifelse(TMB >= 2*870.942857, 'outliers', 
                                   ifelse(TMB <= 441.885714,
                                          'first_quartile','others')))) %>%
  mutate(status = with(., ifelse(Tumor_Sample_Barcode %in% smc56_snv$Tumor_Sample_Barcode,
                                 'SMC5/6 SmV', 'No SMC5/6 SmV')))
# percentage of tumors with SMC5/6 SmV ####
df_plot %>% group_by(quartile, status) %>%
  summarise(n=n())
# low_quartile = 100*19/(115+19) = 14.1791%
# others = 100*95(95+254) = 27.2206%
# outliers = 100*57/(5+57) = 91.9355%
# show TMB of all Hartwig patients 
# show the quartile in different colors
df_plot %>%
  ggplot(aes(x=Tumor_Sample_Barcode,y=TMB,color=quartile)) +
  geom_point() +
  geom_vline(xintercept = c('WIDE01010921T','DRUP01070146T'),
             linetype = 'longdash') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = 'top') +
  scale_color_manual(values = c('#5773CCFF','#5773CCFF','#FFB900FF'))
ggsave(paste(output_dir, '/figures/all_pass_tmb.eps', sep = ''), dpi = 320)
# highlight samples with SMC5/6 SmV
df_plot %>%
  ggplot(aes(x=Tumor_Sample_Barcode,y=TMB,color=status)) +
  geom_point() +
  geom_vline(xintercept = c('WIDE01010921T','DRUP01070146T'),
             linetype = 'longdash') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = 'top') +
  scale_color_manual(values = c('black', 'red'))
ggsave(paste(output_dir, '/figures/highlight_all_pass_tmb.eps', sep = ''), dpi = 320)
# different type of SMC5/6 mutations ----
cohort_master <- list.files(paste(output_dir, '/data', sep = ''), 
                            full.names = T, pattern = 'cohort.csv$') %>%
  map(read.csv, header=T) %>%
  bind_rows()
tmb_master <- read.csv(paste(output_dir, '/data/pass_tmb_master.csv', sep = ''), header = T)
smc56_smv_tumors <- filter(cohort_master, Hugo_Symbol %in% smc56_gene & 
                             alteration_type %in% c('del_mutation','nsnd_mutation','syn_mutation'))
## based on alteration_type, not COHORT ----
df_plot <- tmb_master %>% 
  filter(Tumor_Sample_Barcode %in% condense_survival$sampleId) %>% # show only those with complete clinical info
  left_join(unique(smc56_smv_tumors[,c('Tumor_Sample_Barcode','alteration_type')]), 
             by=join_by(Tumor_Sample_Barcode), relationship='one-to-many') %>%
  mutate(Tumor_Sample_Barcode = fct_reorder(Tumor_Sample_Barcode, -desc(TMB))) %>%
  mutate_at(c('alteration_type'), replace_na, 'others') %>%
  mutate(status = with(., ifelse(alteration_type=='del_mutation','SMC5/6 del',
                                 ifelse(alteration_type=='nsnd_mutation','SMC5/6 nsnd',
                                        ifelse(alteration_type=='syn_mutation','SMC5/6 syn','others')))))
df_plot %>%
  mutate(Tumor_Sample_Barcode = fct_reorder(Tumor_Sample_Barcode, -desc(TMB))) %>%
  ggplot(aes(x=Tumor_Sample_Barcode,y=TMB,color=status)) +
  geom_point() +
  geom_vline(xintercept = c('WIDE01010921T','DRUP01070146T'),
             linetype = 'longdash') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = 'top') +
  scale_color_manual(values = c('SMC5/6 del'='#A52A2A', 'SMC5/6 nsnd'='#80B19B', 'others'='#CCCCCC', 'SMC5/6 syn'='#FAA0A0'))
ggsave(paste(output_dir, '/figures/alteration_type_all_tmb.eps', sep = ''), 
       dpi = 320)
# calculate the percentage of tumors that carry different types of
# SMC5/6 mutations in each of the highlighted section
tmb_master %>%
  mutate(quartile=with(.,ifelse(TMB>=2*10.88678571,'outliers',
                                ifelse(TMB<=5.52357143,'first_quartile','others')))) %>%
  group_by(quartile) %>%
  summarise(n=n())

unique(df_plot[,c('Tumor_Sample_Barcode','TMB','alteration_type')]) %>% 
  mutate(quartile=with(.,ifelse(TMB>=2*10.88678571,'outliers',
                                ifelse(TMB<=5.52357143,'first_quartile','others')))) %>%
  group_by(quartile,alteration_type) %>%
  summarise(n=n())

# quartile  n
# first_quartile  182
# others  468
# outliers  76

# quartile
## based on COHORT, not alteration_type ----
df_plot<- tmb_master %>% 
  filter(Tumor_Sample_Barcode %in% condense_survival$sampleId) %>% # show only those with complete clinical info
  inner_join(unique(cohort_master[,c('Tumor_Sample_Barcode','COHORT')]), 
             by=join_by(Tumor_Sample_Barcode), relationship='one-to-many') %>%
  mutate(Tumor_Sample_Barcode = fct_reorder(Tumor_Sample_Barcode, -desc(TMB))) %>%
  mutate(status = with(., ifelse(COHORT=='del','SMC5/6 del',
                                 ifelse(COHORT=='nsnd','SMC5/6 nsnd',
                                        ifelse(COHORT=='syn','SMC5/6 syn','others')))))
df_plot %>%
  mutate(Tumor_Sample_Barcode = fct_reorder(Tumor_Sample_Barcode, -desc(TMB))) %>%
  ggplot(aes(x=Tumor_Sample_Barcode,y=TMB,color=status)) +
  geom_point() +
  geom_vline(xintercept = c('WIDE01010921T','DRUP01070146T'),
             linetype = 'longdash') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = 'top') +
  scale_color_manual(values = c('SMC5/6 del'='#A52A2A', 'SMC5/6 nsnd'='#80B19B', 'others'='#CCCCCC', 'SMC5/6 syn'='#FAA0A0'))
ggsave(paste(output_dir, '/figures/cohort_all_tmb.eps', sep = ''), 
       dpi = 320)
