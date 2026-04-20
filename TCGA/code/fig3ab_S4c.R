# Fig3a, 3b, and FigS4c
# ENVIRONMENT ####

library('ggplot2')
library('tidyverse')
library('ggpubr')

input_dir <- 'TCGA'
output_dir <- 'TCGA'
# INPUT ####
gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
          'NSMCE4A', 'NSMCE4B', 'EID3',
          'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')
group_color <- c('#FAA0A0', '#A52A2A', '#80B19B', '#C6B3DD', '#702963', '#28446f')
names(group_color) <- c('syn', 'del', 'nonsyn_nondel','low_cnv', 'high_cnv', 'wt')

#i believe the clinical data have aneuploidy score stored 
survival_files <- list.files(path = dir('/cBioPortal/pan_can_atlas',
                                        full.names = TRUE),
                             full.names = TRUE, pattern = 'data_clinical_sample.txt')
survival_ploidy <- purrr::map(survival_files, read.delim, 
                              sep = '\t', header = TRUE, comment.char = '#') %>%
  bind_rows()
cohort_sample <- list.files(path = input_dir, 
  full.names = TRUE, pattern = 'cohort_for_fpkm.csv$') %>%
  purrr::map(read.csv, header = TRUE)
colnames(cohort_sample[[2]]) <- colnames(cohort_sample[[4]])
#keep the primary tumors only 
master_cohort_sample <- bind_rows(cohort_sample) %>%
  filter(stringr::str_detect(Tumor_Sample_Barcode, 'A$') == TRUE) %>%
  mutate(Tumor_Sample_Barcode = stringr::str_extract(Tumor_Sample_Barcode, '^.{15}')) %>%
  select(any_of(c('COHORT', 'Tumor_Sample_Barcode', 'project'))) %>%
  unique()
df <- inner_join(master_cohort_sample, survival_ploidy, by = join_by(Tumor_Sample_Barcode == SAMPLE_ID)) %>%
  filter(COHORT %in% c('syn', 'del', 'nonsyn_nondel', 'low_cnv', 'high_cnv', 'wt')) %>%
  unique()
df$COHORT <- factor(df$COHORT, levels = c('syn', 'del', 'nonsyn_nondel', 'low_cnv', 'high_cnv', 'wt'))
# S4c ----
df %>%
  ggplot(aes(x = COHORT, y = ANEUPLOIDY_SCORE, color = COHORT)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'bottom',
    legend.title = element_blank()
  ) +
  scale_color_manual(values = group_color) +
  labs(x = NULL, y = 'Aneuploidy score') + 
  stat_compare_means(comparisons = list(c('syn', 'wt'), c('del', 'wt'), c('nonsyn_nondel', 'wt'),
                                       c('low_cnv', 'wt'), c('high_cnv', 'wt')),
                     method = 'wilcox.test', hide.ns = TRUE, label = 'p')
ggsave(paste(output_dir, '/figures/abs_cnv_aneuploidy.eps', sep = ''),
       dpi = 320, width = 5.5, height = 4.5, units = 'in')
pairwise.wilcox.test(df$ANEUPLOIDY_SCORE, df$COHORT)


## show the SNV and CNV groups separately 
#load the input file
cohort_sample <- 
  list.files(path = input_dir, full.names = TRUE, pattern = 'cohort_for_fpkm.csv$') %>%
  purrr::map(read.csv, header = TRUE) 
colnames(cohort_sample[[2]]) <- colnames(cohort_sample[[4]])
cohort_sample <- bind_rows(cohort_sample)

#show the SNV
df %>%
  filter(COHORT %in% c('syn', 'del', 'nonsyn_nondel', 'wt')) %>%
  ggplot(aes(x = COHORT, y = ANEUPLOIDY_SCORE, fill = COHORT)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'bottom',
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = c('#FAA0A0', '#A52A2A', '#80B19B', '#28446f')) +
  labs(x = NULL, y = 'Aneuploidy score')
pairwise.wilcox.test(df$ANEUPLOIDY_SCORE, df$COHORT, p.adjust.method = 'none')
ggsave(paste(output_dir, '/figures/separated_snv_aneuploidy_abs_cnv.eps',
             sep = ''), dpi = 320)

#show the CNV
df %>%
  filter(COHORT %in% c('low_cnv', 'high_cnv', 'wt')) %>%
  ggplot(aes(x = COHORT, y = ANEUPLOIDY_SCORE, fill = COHORT)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'bottom',
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = group_color) +
  labs(x = NULL, y = 'Aneuploidy score')
pairwise.wilcox.test(df$ANEUPLOIDY_SCORE, df$COHORT, p.adjust.method = 'none')
ggsave(paste(output_dir, '/figures/separated_cnv_aneuploidy_abs_cnv.eps',
             sep = ''), dpi = 320)
