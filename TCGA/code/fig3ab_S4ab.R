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
group_color <- c('#FAA0A0', '#A52A2A', '#80B19B', '#C6B3DD', '#702963', '#28446f')
names(group_color) <- c('syn', 'del', 'nonsyn_nondel','low_cnv', 'high_cnv', 'wt')
## ANALYSIS ----

### mutational burden ----
tmb_master <- read.csv(paste(input_dir, '/gdc_tmb.csv', sep = ''), header = TRUE)
df <- inner_join(tmb_master, cohort_master, by=join_by(Tumor_Sample_Barcode)) %>%
  select(all_of(c('Tumor_Sample_Barcode','COHORT','TMB'))) %>%
  unique()
### 3a - altered vs non-altered ----
df %>% 
  mutate(color=with(., ifelse(COHORT == 'wt','non-altered', 'altered'))) %>%
  select(all_of(c('color','TMB','Tumor_Sample_Barcode'))) %>%
  unique() %>%
  ggplot(aes(x=color, y=log10(TMB))) +
  geom_boxplot() +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = 'bottom',
    legend.title = element_blank()
  ) +
  stat_compare_means(method = 'wilcox', label = 'p.signif')
### 3b - SmV, CNV, and non-altered ----
df %>% 
  mutate(color=with(., ifelse(COHORT=='wt','non-altered',
                              ifelse(COHORT %in% c('low_cnv','high_cnv'),'CNV','SmV')))) %>%
  select(all_of(c('color','Tumor_Sample_Barcode','TMB'))) %>%
  unique() %>%
  mutate(color=fct_relevel(color, 'SmV','CNV','non-altered')) %>%
  ggplot(aes(x=color,y=log10(TMB))) +
  geom_boxplot() +
  theme_bw() + 
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = 'bottom',
    legend.title = element_blank()
  ) +
  stat_compare_means(comparisons = list(c('SmV','CNV'),c('CNV','non-altered'),c('SmV','non-altered')), 
                     method = 'wilcox', label = 'p.signif')

### S4a - types of CNV ----
df %>% 
  filter(COHORT %in% c('low_cnv', 'high_cnv', 'wt')) %>%
  mutate(COHORT = forcats::fct_relevel(COHORT, 'low_cnv', 'high_cnv', 'wt')) %>%
  ggplot(aes(x = COHORT, y = log10(TMB), color = COHORT)) +
  #geom_pwc(method = 'wilcox.test', hide.ns = 'p') +
  geom_boxplot() +
  theme_bw() + 
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = 'bottom',
    legend.title = element_blank()
  ) +
  scale_color_manual(values = group_color) +
  stat_compare_means(comparisons = list(c('low_cnv','wt'),c('high_cnv','wt')),
                     method = 'wilcox', label = 'p.signif')
ggsave(paste(output_dir, '/figures/abs_cnv_tmb_cnv.eps', sep = ''),
       dpi = 320)
### S4b - show them as overlapping region ----
df %>% 
  filter(COHORT %in% c('low_cnv', 'high_cnv', 'wt')) %>%
  mutate(COHORT = fct_relevel(COHORT,'wt','high_cnv','low_cnv')) %>%
  ggplot(aes(x = log10(TMB), y = COHORT, fill = COHORT)) +
  #geom_pwc(method = 'wilcox.test', hide.ns = 'p') +
  geom_density_ridges(scale = 1.25) +
  theme_ridges() + 
  theme(
    legend.position = 'none'
  ) +
  scale_fill_manual(values = group_color) 
ggsave(paste(output_dir, '/figures/density_abs_cnv_tmb_cnv.eps', sep = ''), dpi = 320)
