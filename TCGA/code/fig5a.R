# ENVIRONMENT ----
library('tidyverse')
library('ggplot2')
library('MutationalPatterns')
library('ggpubr')
input_dir <- 'TCGA'
output_dir <- 'TCGA'
color_code <- c('#DAF7A6','#FFC300',
                '#EEF3FFFF','#CDE5F9FF',
                '#85B6CEFF',
                '#64A8A8FF','#4A9152FF','#45681EFF','#54450FFF','#362904FF',
                c(rep('#dddddc',50)))
names(color_code) <- c('SBS2','SBS13', #apobec
                       'SBS10a', 'SBS10b', #PolE mutation
                       'SBS14', #PolE+MMRd
                       'SBS6','SBS15','SBS21','SBS26','SBS44', #MMRd
                       'SBS84','SBS85', #AID
                       'SBS1','SBS5', #Clock-like
                       'SBS4','SBS29','SBS92', #tobaco smoking
                       'SBS11','SBS25','SBS31','SBS35','SBS86','SBS87',#chemotherapy
                       "SBS3","SBS7a","SBS7b","SBS7c","SBS7d","SBS8",
                       "SBS9","SBS10c","SBS10d","SBS12","SBS16","SBS17a",
                       "SBS17b","SBS18","SBS19","SBS20","SBS22","SBS23",
                       "SBS24","SBS28","SBS30","SBS32","SBS33","SBS34",
                       "SBS36","SBS37","SBS38","SBS39","SBS40","SBS41",
                       "SBS42","SBS88","SBS89","SBS90","SBS91","SBS93","SBS94")

## INPUT ----
gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
          'NSMCE4A', 'NSMCE4B', 'EID3',
          'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')
### list of point mutations ----
snv_master <- read.csv(paste(input_dir, '/data/gdc_snv.csv', sep = ''), header = TRUE)
### 96-trinucleotide context ----
exported_trinu <- read.csv(paste(input_dir,'/data/separated_master_96.csv',sep=''),
                           header=T)
### overlapping cohorts ----
cohort_sample <- 
  list.files(path = input_dir, full.names = TRUE, pattern = 'cohort_for_fpkm.csv$') %>%
  purrr::map(read.csv, header = TRUE) 
colnames(cohort_sample[[2]]) <- colnames(cohort_sample[[4]])
cohort_master <- bind_rows(cohort_sample) 
# ANALYSIS ----
del_tumors <- filter(cohort_master, COHORT == 'del')
syn_tumors <- filter(cohort_master, COHORT == 'syn')
nonsyn_nondel_tumors <- filter(cohort_master, COHORT == 'nonsyn_nondel')
low_tumors <- filter(cohort_master, COHORT == 'low_cnv')
high_tumors <- filter(cohort_master, COHORT == 'high_cnv')
wt_tumors <- filter(cohort_master, COHORT == 'wt')
del_spectra <- select(ungroup(exported_trinu), all_of(c(del_tumors$Tumor_Sample_Barcode)))
syn_spectra <- select(ungroup(exported_trinu), all_of(c(syn_tumors$Tumor_Sample_Barcode)))
nonsyn_nondel_spectra <- select(ungroup(exported_trinu), all_of(c(nonsyn_nondel_tumors$Tumor_Sample_Barcode)))
low_spectra <- select(ungroup(exported_trinu), any_of(c(low_tumors$Tumor_Sample_Barcode)))
high_spectra <- select(ungroup(exported_trinu), any_of(c(high_tumors$Tumor_Sample_Barcode)))
wt_spectra <- select(ungroup(exported_trinu), any_of(c(wt_tumors$Tumor_Sample_Barcode)))

df <- lapply(list(syn_spectra, del_spectra, nonsyn_nondel_spectra,
                  low_spectra, high_spectra, wt_spectra),
             rowSums) %>%
  lapply(., as.data.frame)
colnames(df[[1]]) <- 'syn'
colnames(df[[2]]) <- 'del'
colnames(df[[3]]) <- 'nonsyn_nondel'
colnames(df[[4]]) <- 'low_cnv'
colnames(df[[5]]) <- 'high_cnv'
colnames(df[[6]]) <- 'wt'
df <- cbind(df[[1]], df[[2]], df[[3]], df[[4]], df[[5]], df[[6]])

### fit to signature -----
cosmic_v3.2 <- get_known_signatures(source = 'COSMIC_v3.2')
df_fitting <- fit_to_signatures(df, cosmic_v3.2)
df_fitting$contribution %>% as.data.frame() %>%
  write.csv(., paste(input_dir, '/data/sbs_sig.csv', sep = ''))

df_rel <- apply(df_fitting$contribution, 2, function(x) x / sum(x)) %>%
  as_tibble(rownames = NA)
df_plot <- df_rel %>%
  select(all_of(c('syn','del','nonsyn_nondel','wt'))) %>%
  mutate('sbs' = rownames(df)) %>%
  pivot_longer(!c(sbs), names_to = 'group', values_to = 'rel_cont') %>%
  mutate(group = forcats::fct_relevel(group, 'del', 'syn', 'nonsyn_nondel','wt')) %>%
  mutate(sbs = forcats::fct_relevel(sbs, 'SBS2','SBS13', #apobec
                                    'SBS10a', 'SBS10b', #PolE mutation
                                    'SBS14', #PolE+MMRd
                                    'SBS6','SBS15','SBS21','SBS26','SBS44', #MMRd
                                    'SBS84','SBS85', #AID
                                    'SBS1','SBS5', #Clock-like
                                    'SBS4','SBS29','SBS92', #tobaco smoking
                                    'SBS11','SBS25','SBS31','SBS35','SBS86','SBS87',#chemotherapy
                                    "SBS3","SBS7a","SBS7b","SBS7c","SBS7d","SBS8",
                                    "SBS9","SBS10c","SBS10d","SBS12","SBS16","SBS17a",
                                    "SBS17b","SBS18","SBS19","SBS20","SBS22","SBS23",
                                    "SBS24","SBS28","SBS30","SBS32","SBS33","SBS34",
                                    "SBS36","SBS37","SBS38","SBS39","SBS40","SBS41",
                                    "SBS42","SBS88","SBS89","SBS90","SBS91","SBS93","SBS94")) 
ggplot(df_plot, aes(x=group, y=rel_cont, fill=sbs)) +
  geom_col(position='stack') + 
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  scale_fill_manual(values = color_code) +
  labs(x='SMC5/6 gene status', y='Relative contribution')
ggsave(paste(output_dir, '/figures/SBSsig_unified_4groups.eps', sep = ''),
       dpi = 320, width = 7, height = 5, units = 'in')

