# ENVIRONMENT ####
library('tidyverse')
library('ggplot2')
library('UpSetR')
working_dir <- 'PCAWG'
# INPUT ####
smc56_snv <- read.csv(paste(working_dir, '/data/smc56_snv.csv', sep = ''),
                      header = TRUE)
smc56_low_cnv <- read.csv(paste(working_dir, '/data/smc56_low_cnv.csv', sep = ''),
                          header = TRUE)
smc56_high_cnv <- read.csv(paste(working_dir, '/data/smc56_high_cnv.csv', sep = ''),
                           header = TRUE)
# to keep the whitelist patients only
white_list_sample <- read.csv(paste(working_dir, '/data/white_list_sample_sheet.csv', sep = ''),
                              header = TRUE)
# to keep the primary samples only 
clinical_sample <- read.delim('cBioPortal/pancan_pcawg_2020/data_clinical_sample.txt',
                              header = TRUE, sep = '\t', comment.char = '#', fill = FALSE, quote = '')
# ANALYSIS ####
primary_sample <- filter(clinical_sample, SAMPLE_TYPE=='Primary') %>%
  inner_join(white_list_sample, by = join_by(SAMPLE_ID==icgc_specimen_id))

# make a binary data table that is compatible to the upsetR command ####
x2 <- select(smc56_low_cnv, all_of(c('Gene_Symbol', 'Tumor_Sample_Barcode'))) %>%
  filter(Tumor_Sample_Barcode %in% primary_sample$aliquot_id) %>%
  unique() %>%
  mutate(alteration_type = 'Low CNV')
x3 <- select(smc56_high_cnv, all_of(c('Gene_Symbol', 'Tumor_Sample_Barcode'))) %>%
  filter(Tumor_Sample_Barcode %in% primary_sample$aliquot_id) %>%
  unique() %>%
  mutate(alteration_type = 'High CNV')

## remove the intron mutations from the analysis ----
x1 <- filter(smc56_snv, Variant_Classification != 'Intron') %>%
  select(all_of(c('Hugo_Symbol', 'Tumor_Sample_Barcode'))) %>%
  filter(Tumor_Sample_Barcode %in% primary_sample$aliquot_id) %>%
  unique() %>%
  mutate(alteration_type = 'SmV')
colnames(x1) <- c('Gene_Symbol', 'Tumor_Sample_Barcode', 'alteration_type')
# each tumor + SMC56 gene was counted only once
df <- rbind(x1[,c('Gene_Symbol', 'Tumor_Sample_Barcode')], 
            x2[,c('Gene_Symbol', 'Tumor_Sample_Barcode')],
            x3[,c('Gene_Symbol', 'Tumor_Sample_Barcode')]) %>%
  unique() %>%
  as.data.frame() %>%
  mutate(Gene_Symbol = with(., ifelse(Gene_Symbol=='ANKRD32',
                                      'SLF1', ifelse(Gene_Symbol=='FAM178A',
                                                     'SLF2', ifelse(Gene_Symbol=='EID3',
                                                                    'NSMCE4B', ifelse(Gene_Symbol=='NDNL2',
                                                                                      'NSMCE3', Gene_Symbol)))))) %>%
  mutate_at(c('Gene_Symbol', 'Tumor_Sample_Barcode'), as.character) %>%
  table() %>%
  as.data.frame() %>%
  pivot_wider(names_from = Gene_Symbol, values_from = Freq) %>%
  as.data.frame()
# plot the UpSet plot
upset(df, sets = c('SLF2','SLF1','NSMCE4B','NSMCE4A','NSMCE3',
                   'NSMCE2','NSMCE1','SMC6','SMC5'), keep.order = TRUE,
      mainbar.y.label = "Number of tumors", sets.x.label = "Number of tumors",
      order.by = 'freq', nintersects = 50, mb.ratio = c(0.6, 0.4),
      shade.alpha = 1, matrix.dot.alpha = 1
)
