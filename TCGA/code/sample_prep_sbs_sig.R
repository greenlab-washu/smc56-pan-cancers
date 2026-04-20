library('tidyverse')
library('ggplot2')
library('MutationalPatterns')
library('ggpubr')
library('ggridges')

input_dir <- 'TCGA'
output_dir <- 'TCGA'

## INPUT ----
gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
          'NSMCE4A', 'NSMCE4B', 'EID3',
          'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')
### list of point mutations ----
snv_master <- read.csv(paste(input_dir, '/gdc_snv.csv', sep = ''), header = TRUE)
### overlapping cohorts ----
cohort_sample <- 
  list.files(path = input_dir, full.names = TRUE, pattern = 'cohort_for_fpkm.csv$') %>%
  purrr::map(read.csv, header = TRUE) 
colnames(cohort_sample[[2]]) <- colnames(cohort_sample[[4]])
cohort_master <- bind_rows(cohort_sample) 

# trinucleotide context ----
trinu <- select(snv_master, all_of(c('Tumor_Sample_Barcode', 'CONTEXT', 'HGVSc'))) %>% 
  mutate('Mutation_type' = stringr::str_extract(HGVSc, '[:upper:]>[:upper:]')) %>% 
  mutate('Trinucleotide' = stringr::str_sub(CONTEXT, 5, 7)) %>%
  filter(!is.na(Mutation_type)) 
trinu$Mutation_type <- stringr::str_replace_all(trinu$Mutation_type, c('A>C' = 'T>G', 'A>T' = 'T>A', 'A>G' = 'T>C',
                                                                       'G>C' = 'C>G', 'G>T' = 'C>A', 'G>A' = 'C>T'))
trinu$Trinucleotide <- stringr::str_replace_all(trinu$Trinucleotide, c('TGT'="ACA", 'TGG'="ACC", 'TGC'="ACG", 'TGA'="ACT",
                                                                       'GGT'="CCA", 'GGG'="CCC", 'GGC'="CCG", 'GGA'="CCT",
                                                                       'CGT'="GCA", 'CGG'="GCC", 'CGC'="GCG", 'CGA'="GCT",
                                                                       'AGT'="TCA", 'AGG'="TCC", 'AGC'="TCG", 'AGA'="TCT",
                                                                       'TAT'="ATA", 'TAG'="ATC", 'TAC'="ATG", 'TAA'="ATT",
                                                                       'GAT'="CTA", 'GAG'="CTC", 'GAC'="CTG", 'GAA'="CTT",
                                                                       'CAT'="GTA", 'CAG'="GTC", 'CAC'="GTG", 'CAA'="GTT",
                                                                       'AAT'="TTA", 'AAG'="TTC", 'AAC'="TTG", 'AAA'="TTT"))
exported_trinu <- trinu %>%
  group_by(Mutation_type, Trinucleotide, Tumor_Sample_Barcode) %>%
  summarise(n = n()) %>% 
  pivot_wider(names_from = Tumor_Sample_Barcode, values_from = n) %>%
  replace(is.na(.), 0)
exported_trinu <- exported_trinu[order(exported_trinu$Mutation_type, exported_trinu$Trinucleotide),]
write.csv(exported_trinu, paste(output_dir, '/data/separated_master_96.csv', sep = ''),
          row.names = FALSE)

