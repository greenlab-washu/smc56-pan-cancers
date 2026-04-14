# Oct 23rd, 2025
# the current cohort files are too confusing, so I will reorganize them 
# also to reflect the fact that a tumor belongs to one cohort
# can also carry other different types of alterations 

# ENVIRONMENT ----
library('tidyverse')
input_dir <- '/storage2/fs1/abby.green/Active/Hartwig'
output_dir <- '/storage2/fs1/abby.green/Active/Users/thi/SMC56_de_novo/2025/Hartwig'
smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
                'NSMCE4A', 'NSMCE4B', 'EID3',
                'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')
# INPUT ----
pass_smv_master <- list.files(paste(output_dir, '/data', sep = ''), 
                              full.names = T, pattern = 'patients_pass_only') %>% 
  map(read.csv, header=T) %>% 
  bind_rows()
cnv_master <- read.csv(paste(output_dir, '/data/all_patients_cnv_master.csv', 
                            sep = ''), header=T)
# defining SmV cohorts ----
# based on the non-intronic mutations only 
non_intron_smv_master <- filter(pass_smv_master, Variant_Classification != 'Intron') %>%
  mutate(alteration_type=with(., ifelse(Variant_Classification=='Silent',
                                        'syn_mutation', ifelse(IMPACT=='HIGH' | (SIFT!='' & str_extract(SIFT,'(?<=\\().+(?=\\))')<=0.05),
                                                               'del_mutation','nsnd_mutation'))))
del_tumors <- filter(non_intron_smv_master, Hugo_Symbol %in% smc56_gene & alteration_type=='del_mutation') 
nsnd_tumors <- filter(non_intron_smv_master, Hugo_Symbol %in% smc56_gene & alteration_type=='nsnd_mutation') %>%
  filter(!Tumor_Sample_Barcode %in% del_tumors$Tumor_Sample_Barcode) 
syn_tumors <- filter(non_intron_smv_master, Hugo_Symbol %in% smc56_gene & alteration_type=='syn_mutation') %>%
  filter(!Tumor_Sample_Barcode %in% union(del_tumors$Tumor_Sample_Barcode, nsnd_tumors$Tumor_Sample_Barcode))
# defining CNV cohorts ----
df_cnv <- cnv_master %>%
  mutate(alteration_type=with(., ifelse(maxCopyNumber <= 1.5,
                                        'high_cnv', ifelse(minCopyNumber >= 2.5,
                                                           'low_cnv', 'normal_cnv')))) %>%
  select(all_of(c('Tumor_Sample_Barcode','gene','alteration_type'))) %>%
  unique()
colnames(df_cnv) <- c('Tumor_Sample_Barcode','Hugo_Symbol','alteration_type')
low_cnv_tumors <- filter(df_cnv, Hugo_Symbol %in% smc56_gene & alteration_type=='low_cnv')
high_cnv_tumors <- filter(df_cnv, Hugo_Symbol %in% smc56_gene & alteration_type=='high_cnv')
wt_tumors <- filter(df_cnv, 
                    !Tumor_Sample_Barcode %in% union(del_tumors$Tumor_Sample_Barcode,union(
                      nsnd_tumors$Tumor_Sample_Barcode,union(
                        syn_tumors$Tumor_Sample_Barcode, union(
                          low_cnv_tumors$Tumor_Sample_Barcode, high_cnv_tumors$Tumor_Sample_Barcode))))) %>%
  mutate(COHORT='wt')
del_cnv_tumors <- filter(df_cnv, Tumor_Sample_Barcode %in% del_tumors$Tumor_Sample_Barcode)
nsnd_cnv_tumors <- filter(df_cnv, Tumor_Sample_Barcode %in% nsnd_tumors$Tumor_Sample_Barcode)
syn_cnv_tumors <- filter(df_cnv, Tumor_Sample_Barcode %in% syn_tumors$Tumor_Sample_Barcode)
# output ----
wt_cohort <- filter(non_intron_smv_master, Tumor_Sample_Barcode %in% wt_tumors$Tumor_Sample_Barcode) %>%
  select(all_of(c('Tumor_Sample_Barcode','Hugo_Symbol','alteration_type'))) %>%
  unique() %>%
  mutate(COHORT='wt') %>%
  rbind(wt_tumors)
del_cohort <- filter(non_intron_smv_master, Tumor_Sample_Barcode %in% del_tumors$Tumor_Sample_Barcode) %>%
  select(all_of(c('Tumor_Sample_Barcode','Hugo_Symbol','alteration_type'))) %>%
  unique() %>%
  bind_rows(del_cnv_tumors) %>%
  mutate(COHORT='del') 
nsnd_cohort <- filter(non_intron_smv_master, Tumor_Sample_Barcode %in% nsnd_tumors$Tumor_Sample_Barcode) %>%
  select(all_of(c('Tumor_Sample_Barcode','Hugo_Symbol','alteration_type'))) %>%
  unique() %>%
  bind_rows(nsnd_cnv_tumors) %>%
  mutate(COHORT='nsnd')
syn_cohort <- filter(non_intron_smv_master, Tumor_Sample_Barcode %in% syn_tumors$Tumor_Sample_Barcode) %>%
  select(all_of(c('Tumor_Sample_Barcode','Hugo_Symbol','alteration_type'))) %>%
  unique() %>%
  bind_rows(syn_cnv_tumors) %>%
  mutate(COHORT='syn')
# Overlapping tumors between low vs high CNV groups will be omitted 
# to keep it in consistent with the TCGA analysis
low_cnv_cohort <- filter(non_intron_smv_master,
                         Tumor_Sample_Barcode %in% setdiff(low_cnv_tumors$Tumor_Sample_Barcode, high_cnv_tumors$Tumor_Sample_Barcode)) %>%
  select(all_of(c('Tumor_Sample_Barcode','Hugo_Symbol','alteration_type'))) %>%
  unique() %>%
  bind_rows(low_cnv_tumors) %>%
  mutate(COHORT='low_cnv') %>%
  filter(Tumor_Sample_Barcode %in% setdiff(low_cnv_tumors$Tumor_Sample_Barcode, high_cnv_tumors$Tumor_Sample_Barcode))
high_cnv_cohort <- filter(non_intron_smv_master,
                          Tumor_Sample_Barcode %in% setdiff(high_cnv_tumors$Tumor_Sample_Barcode, low_cnv_tumors$Tumor_Sample_Barcode)) %>%
  select(all_of(c('Tumor_Sample_Barcode','Hugo_Symbol','alteration_type'))) %>%
  unique() %>%
  bind_rows(high_cnv_tumors) %>%
  mutate(COHORT='high_cnv') %>%
  filter(Tumor_Sample_Barcode %in% setdiff(high_cnv_tumors$Tumor_Sample_Barcode, low_cnv_tumors$Tumor_Sample_Barcode))

write_csv(wt_cohort, paste(output_dir, '/data/smc56_wt_cohort.csv', sep = ''))
write_csv(syn_cohort, paste(output_dir, '/data/smc56_syn_cohort.csv', sep = ''))
write_csv(del_cohort, paste(output_dir, '/data/smc56_del_cohort.csv', sep = ''))
write_csv(nsnd_cohort, paste(output_dir, '/data/smc56_nsnd_cohort.csv', sep = ''))
write_csv(low_cnv_cohort, paste(output_dir, '/data/smc56_low_cnv_cohort.csv', sep = ''))
write_csv(high_cnv_cohort, paste(output_dir, '/data/smc56_high_cnv_cohort.csv', sep = ''))
