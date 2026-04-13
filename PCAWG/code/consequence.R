# July 1st, 2025  
# Thi Tran 
# consequence of SMC5/6 SmV
# based on types, not by deleterious prediction

# ENVIRONMENT ####
library('tidyverse')
library('ggplot2')
working_dir <- '/storage2/fs1/abby.green/Active/Users/thi/SMC56_de_novo/2025/PCAWG'
theme_general <-   theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_blank(),
    legend.position = 'bottom'
  ) 
# INPUT ####
raw_smc56_smv <- read.csv(paste(working_dir, '/data/smc56_snv.csv', sep = ''),
                      header = TRUE)
# to keep the whitelist patients only
white_list_sample <- read.csv(paste(working_dir, '/data/white_list_sample_sheet.csv', sep = ''),
                              header = TRUE)
# to keep the primary samples only 
clinical_sample <- read.delim('/storage2/fs1/abby.green/Active/RawData/TCGA_PCAWG/cBioPortal/pancan_pcawg_2020/data_clinical_sample.txt',
                              header = TRUE, sep = '\t', comment.char = '#', fill = FALSE, quote = '')
# ANALYSIS ####
primary_sample <- filter(clinical_sample, SAMPLE_TYPE=='Primary') %>%
  inner_join(white_list_sample, by = join_by(SAMPLE_ID==icgc_specimen_id))
smc56_smv <- filter(raw_smc56_smv, Tumor_Sample_Barcode %in% primary_sample$aliquot_id) 
smc56_smv %>% 
  select(all_of(c('Hugo_Symbol','Tumor_Sample_Barcode','Variant_Classification'))) %>%
  unique() %>%
  mutate(gene_name=with(., ifelse(Hugo_Symbol=='ANKRD32',
                                  'SLF1', ifelse(Hugo_Symbol=='FAM178A',
                                                 'SLF2', ifelse(Hugo_Symbol=='EID3',
                                                                'NSMCE4B', ifelse(Hugo_Symbol=='NDNL2',
                                                                                  'NSMCE3', Hugo_Symbol)))))) %>%
  mutate(gene_name=forcats::fct_relevel(gene_name, 'SMC5','SMC6','NSMCE1',
                                        'NSMCE2','NSMCE3','NSMCE4A',
                                        'NSMCE4B','SLF1','SLF2')) %>%
  ggplot(aes(x=gene_name,fill=Variant_Classification)) +
  geom_bar() + 
  theme_general +
  ylab('Number of tumors') +
  xlab('SMC5/6 genes')
ggsave(paste(working_dir, '/figures/consequence.eps', sep = ''),
       dpi = 320, width = 6.8, height = 4, units = 'in')

# omit the intron mutations ----
smc56_smv %>% 
  select(all_of(c('Hugo_Symbol','Tumor_Sample_Barcode','Variant_Classification'))) %>%
  filter(Variant_Classification != 'Intron') %>%
  unique() %>%
  mutate(gene_name=with(., ifelse(Hugo_Symbol=='ANKRD32',
                                  'SLF1', ifelse(Hugo_Symbol=='FAM178A',
                                                 'SLF2', ifelse(Hugo_Symbol=='EID3',
                                                                'NSMCE4B', ifelse(Hugo_Symbol=='NDNL2',
                                                                                  'NSMCE3', Hugo_Symbol)))))) %>%
  mutate(gene_name=forcats::fct_relevel(gene_name, 'SMC5','SMC6','NSMCE1',
                                        'NSMCE2','NSMCE3','NSMCE4A',
                                        'NSMCE4B','SLF1','SLF2')) %>%
  ggplot(aes(x=gene_name,fill=Variant_Classification)) +
  geom_bar() +
  theme_general +
  ylab('Number of tumors') +
  xlab('SMC5/6 genes') +
  scale_fill_brewer(palette = 'PiYG')
ggsave(paste(working_dir, '/figures/no_intron_consequence.eps', sep = ''),
       dpi = 320, width = 6.8, height = 4, units = 'in')

## VAF as well ----
smc56_smv %>% 
  filter(Variant_Classification != 'Intron') %>%
  select(all_of(c('Hugo_Symbol','Tumor_Sample_Barcode','Genome_Change', 'i_VAF'))) %>%
  unique() %>%
  mutate(gene_name=with(., ifelse(Hugo_Symbol=='ANKRD32',
                                  'SLF1', ifelse(Hugo_Symbol=='FAM178A',
                                                 'SLF2', ifelse(Hugo_Symbol=='EID3',
                                                                'NSMCE4B', ifelse(Hugo_Symbol=='NDNL2',
                                                                                  'NSMCE3', Hugo_Symbol)))))) %>%
  mutate(gene_name=forcats::fct_relevel(gene_name, 'SMC5','SMC6','NSMCE1',
                                        'NSMCE2','NSMCE3','NSMCE4A',
                                        'NSMCE4B','SLF1','SLF2')) %>%
  ggplot(aes(x=gene_name,y=as.numeric(i_VAF))) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               binpositions = 'all', dotsize = 0.75) +
  theme_general +
  ylab('Variant allele frequency') +
  ylim(0,1) +
  xlab('SMC5/6 genes')
ggsave(paste(working_dir, '/figures/no_intron_vaf.eps', sep = ''),
       dpi = 320, width = 6.8, height = 4, units = 'in')
