# ENVIRONMENT ----
library('tidyverse')
library('ggplot2')
library('ggpubr')

working_dir <- 'TCGA'

# INPUT ####
## recurrent SMC5/6 mutations ----
smc56_recurrent_master <- read.csv(paste(working_dir, '/data/recurrent_smc56_mutations.csv', sep = ''),
                                row.names = 1)
## load all SNV mutations ----
snv_master <- read.csv(paste(working_dir, '/data/gdc_snv.csv', sep = ''))
## load OVERLAP COHORT ----
overlap_cohort_samples <- list.files(
  path = paste(working_dir, '/data', sep = ''), 
  full.names = TRUE, pattern = 'cohort_for_fpkm.csv$') %>%
  purrr::map(read.csv, header = TRUE)
colnames(overlap_cohort_samples[[2]]) <- colnames(overlap_cohort_samples[[4]])
overlap_cohort_master <- bind_rows(overlap_cohort_samples) %>%
  mutate(patient_id = stringr::str_extract(Tumor_Sample_Barcode, '^.{12}')) %>%
  select(any_of(c('COHORT', 'Tumor_Sample_Barcode', 'project', 'patient_id'))) %>%
  unique()

cohort_sample <- 
  list.files(path = paste(working_dir, '/data', sep = ''), full.names = TRUE, pattern = 'cohort_for_fpkm.csv$') %>%
  purrr::map(read.csv, header = TRUE)
colnames(cohort_sample[[2]]) <- colnames(cohort_sample[[4]])
cohort_sample <- bind_rows(cohort_sample) %>% 
  select(all_of(c('COHORT', 'Tumor_Sample_Barcode')))

tmb_master <- read.csv(paste(input_dir, '/data/gdc_tmb.csv',  sep=''), header = T)
# ANALYSIS ----
df <- inner_join(tmb_master, smc56_recurrent_master,
                 by = join_by(Hugo_Symbol, HGVSp_Short)) %>%
  filter(mutation_type == 'del_mutation' & n > 2) # keep only del and recurrent SMC5/6 mutations
## 2c ----
df %>% 
  ungroup() %>%
  mutate(Hugo_Symbol = forcats::fct_relevel(Hugo_Symbol, 'SMC5', 'SMC6', 'NSMCE4A',
                                            'EID3', 'SLF1')) %>% 
  mutate(HGVSp_Short = forcats::fct_relevel(HGVSp_Short, 'p.R972*', 'p.R978*', 'p.I410Yfs*4', 'p.R279*',
                                            'p.S87L', 'p.R289C', 'p.D89N', 'p.R185Q', 'p.T935A')) %>%
  ggplot(aes(x=HGVSp_Short, y=log10(TMB), color=Hugo_Symbol)) +
  geom_boxplot() +
  theme_bw() + 
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.title = element_blank()
  ) +
  geom_hline(aes(yintercept = median(log10(TMB))), del_tmb, color = '#A52A2A') +
  geom_hline(aes(yintercept = median(log10(TMB))), wt_tmb, color = '#28446f') + 
  scale_x_discrete(labels = c('p.R972* (n=6)', 'p.R978* (n=4)', 'p.I410Yfs*4 (n=30)', 'p.R279* (n=3)',
                              'p.S87L (n=2)', 'p.R289C (n=5)', 'p.D89N (n=3)', 'p.R185Q (n=3)', 'p.T935A (n=3)'))
ggsave(paste(working_dir, '/figures/tmb_del_smc56_hotspot_mut.eps', sep=''),
       width = 7, height = 3.5, dpi = 320, units = 'in')
