# ENVIRONMENT ####
working_dir <- 'TCGA'
library('ggplot2')
library('tidyverse')
library('ggpubr')
gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
          'NSMCE4A', 'NSMCE4B', 'EID3',
          'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')

# INPUT ####
## TMB 
snv_master <- read.csv(paste(working_dir, '/data/gdc_snv.csv', sep = ''))
tmb_master <- read.csv(paste(working_dir, 'data/gdc_tmb.csv', sep = ''), 
                       header = T)
# SMC5/6 mutations vs TMB ####
### number of SMC5/6 SNVs within individual tumors 
smc56_snv_per_tumor <- filter(snv_master, Hugo_Symbol %in% gene) %>%
  filter(stringr::str_detect(Tumor_Sample_Barcode, 'A$') == TRUE) %>% #keep the primary tumors only, whose Tumor_Sample_Barcode ends with an A
  select(all_of(c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'HGVSc', 'Consequence'))) %>%
  unique() %>%
  group_by(Tumor_Sample_Barcode) %>% 
  count() %>%
  ungroup() 
## PLOT ----
### Accumulated SMC5/6 mutations ----
#### scatter plot with best-fit line ----
df <- full_join(tmb_master, smc56_snv_per_tumor, by = join_by(Tumor_Sample_Barcode)) %>%
  unique()
df$n <- replace_na(df$n, 0)
df %>%
  ggplot(aes(x=n, y=TMB)) +
  geom_point() +
  geom_smooth(method = 'lm', color = 'red') +
  stat_regline_equation() +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  stat_regline_equation() +
  labs(x = 'Accumulated number of SMC5/6 mutations', y = 'TMB (number of mutations/Mb)')
ggsave(paste(working_dir, '/figures/smc56_tmb_correlation_scatter.eps', sep = ''),
       dpi = 320, width = 4, height = 3.5, units = 'in')
## Pearson correlation test
cor.test(df$n, df$TMB, alternative = 'greater')
# all statistical information about the best-fit line 
summary(lm(TMB ~ n, data = df))
