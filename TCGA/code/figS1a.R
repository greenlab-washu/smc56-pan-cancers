library('UpSetR')
library('tidyverse')
library('ggplot2')
library('MutationalPatterns')
library('ggpubr')
input_dir <- 'TCGA'
output_dir <- 'TCGA'

# INPUT ####
# Tumor IDs, SMC5/6 genetic alterations
cohort_sample <- 
  list.files(path = paste(input_dir, '/data', sep = ''), full.names = TRUE, pattern = 'cohort_for_fpkm.csv$') %>%
  purrr::map(read.csv, header = TRUE) 
colnames(cohort_sample[[2]]) <- colnames(cohort_sample[[4]])
cohort_master <- bind_rows(cohort_sample)
smc56genes <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3',
                'NSMCE4A', 'EID3', 'SLF1', 'SLF2')
# ANALYSIS ####
df <- filter(cohort_master, COHORT != 'wt' & alteration_type != 'normal_cna' &
               Hugo_Symbol %in% smc56genes & 
               stringr::str_extract(Tumor_Sample_Barcode, '.$') == 'A') 
df_snv <- filter(df, stringr::str_extract(alteration_type, '.{3}$') != 'cna') %>%
  select(all_of(c('Tumor_Sample_Barcode', 'Hugo_Symbol'))) %>%
  unique() %>%
  mutate(freq = 1) %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = freq, values_fill = 0) %>%
  as.data.frame() %>%
  mutate(Tumor_Sample_Barcode = as.character(Tumor_Sample_Barcode))
upset(df_snv, sets = c('SLF2','SLF1','EID3','NSMCE4A','NSMCE3',
                       'NSMCE2','NSMCE1','SMC6','SMC5'), keep.order = TRUE,
      mainbar.y.label = "Number of tumors", sets.x.label = "Number of tumors",
      order.by = 'freq', nintersects = 29, mb.ratio = c(0.6, 0.4),
      shade.alpha = 1, matrix.dot.alpha = 1)