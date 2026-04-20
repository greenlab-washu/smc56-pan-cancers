# ENVIRONMENT ----
library('ggplot2')
library('tidyverse')
library('survminer')
library('survival')

input_dir <- 'TCGA'
output_dir <- 'TCGA'
theme_general <- theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

# INPUT ----
gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
          'NSMCE4A', 'NSMCE4B', 'EID3',
          'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')
overlap_cohort_samples <- 
  list.files(path = input_dir, full.names = TRUE, pattern = 'cohort_for_fpkm.csv$') %>%
  purrr::map(read.csv, header = TRUE) 
colnames(cohort_sample[[2]]) <- c('COHORT', 'alteration_type', 'Tumor_Sample_Barcode',
                                  'Hugo_Symbol', 'project')
overlap_cohort_samples <- bind_rows(overlap_cohort_samples)

snv_master <- read.csv(paste(input_dir, '/gdc_snv.csv', sep = ''), header = TRUE)
## WT cohort ----
wt_cohort_master <- data_frame(COHORT = 'wt', alteration_type = 'wt', 
                               Tumor_Sample_Barcode = unique(overlap_cohort_samples[[6]]$Tumor_Sample_Barcode),
                               Hugo_Symbol = rep('Wild-type', 1984),
                               project = NA
) 

## TMB ----
tmb_master <- select(snv_master, all_of(c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Chromosome', 
                                          'Start_Position', 'End_Position', 'HGVSc'))) %>%
  unique() %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(abs_mut_number = n()) %>%
  mutate(TMB = abs_mut_number/30) %>%
  mutate(logTMB = log10(TMB))
# OVERLAP included 
#different types of genetic alterations
df <- filter(overlap_cohort_master, Hugo_Symbol %in% gene) %>%
  rbind(wt_cohort_master) %>%
  inner_join(tmb_master, by = join_by(Tumor_Sample_Barcode)) %>%
  select(all_of(c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'alteration_type', 'TMB', 'logTMB'))) %>%
  unique() %>%
  mutate(Hugo_Symbol = forcats::fct_relevel(Hugo_Symbol, 'SMC5', 'SMC6', 'NSMCE1',
                                            'NSMCE2', 'NSMCE3', 'NSMCE4A',
                                            'EID3', 'SLF1', 'SLF2', 'Wild-type')) %>%
  mutate(alteration_type = forcats::fct_relevel(alteration_type, 'syn_mutation',
                                                'del_mutation', 'nonsyn_nondel_mutation',
                                                'low_cna', 'high_cna', 'normal_cna', 'wt'))
### SNV ----
#### log10 TMB ----
filter(df, alteration_type %in% c('syn_mutation', 'del_mutation', 'nonsyn_nondel_mutation', 'wt')) %>%
  ggplot(aes(x = Hugo_Symbol, y = logTMB, fill = alteration_type)) +
  geom_boxplot() +
  theme_bw() + 
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = c('#FAA0A0', '#A52A2A', '#80B19B', '#28446f')) +
  labs(x = NULL, y = 'log10(TMB)', fill = NULL)
ggsave(paste(output_dir, '/figures/point_gene_specific_log10_tmb.eps', sep = ''),
       dpi = 320)
### CNV ----
#### log10 TMB ----
filter(df, alteration_type %in% c('low_cna', 'high_cna', 'wt')) %>%
  ggplot(aes(x = Hugo_Symbol, y = logTMB, fill = alteration_type)) +
  geom_boxplot() +
  theme_bw() + 
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = c('#C6B3DD', '#702963', '#28446f')) +
  labs(x = NULL, y = 'log10(TMB)', fill = NULL)
ggsave(paste(output_dir, '/figures/cnv_gene_specific_log10_tmb.eps', sep = ''),
       dpi = 320)