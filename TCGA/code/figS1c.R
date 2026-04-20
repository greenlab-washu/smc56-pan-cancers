# ENVIRONMENT ####
library('ggplot2')
library('tidyverse')

input_dir <- 'TCGA'
output_dir <- 'TCGA'

# INPUT #####
gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
          'NSMCE4A', 'NSMCE4B', 'EID3',
          'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')

groups_color <- c('#B04545', '#609E82')
names(groups_color) <- c('del_mutation', 'syn_mutation')
snv_master <- read.csv(paste(input_dir, '/gdc_snv.csv', sep = ''), header = TRUE) 
#load SNV cohort ID
cohort_sample <- 
  list.files(path = input_dir, full.names = TRUE, pattern = 'cohort_for_fpkm.csv$') %>%
  purrr::map(read.csv, header = TRUE) 
colnames(cohort_sample[[2]]) <- c('COHORT', 'alteration_type', 'Tumor_Sample_Barcode',
                                  'Hugo_Symbol', 'project')
cohort_sample <- bind_rows(cohort_sample)

# ANALYSIS #####
snv_smc <- filter(snv_master, Hugo_Symbol %in% gene) %>%
  mutate_at('Hugo_Symbol', ~ifelse(. == 'NDNL2', 'NSMCE3',
                                   ifelse(. == 'ANKRD32', 'SLF1',
                                          ifelse(. == 'FAM178A', 'SLF2',
                                                 ifelse(. == 'EID3', 'NSMCE4B', .))))) %>%
  mutate(VAF = t_alt_count/(t_alt_count + t_ref_count)) %>%
  mutate(Hugo_Symbol = forcats::fct_relevel(Hugo_Symbol, 'SMC5', 'SMC6', 'NSMCE1', 'NSMCE2',
                                            'NSMCE3', 'NSMCE4A', 'NSMCE4B', 'SLF1', 'SLF2')) %>% 
  mutate(alteration_type = with(., ifelse(stringr::str_detect(HGVSp_Short, '=')==TRUE,
                                          'syn_mutation', ifelse(IMPACT == 'HIGH',
                                                                 'del_mutation', ifelse(SIFT != '' & stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') <= 0.05,
                                                                                        'del_mutation', 'nonsyn_nondel')))))

# dotplot
## all types of small variants
snv_smc %>% 
  ggplot(aes(x=Hugo_Symbol, y=VAF)) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', 
               binpositions = 'all', dotsize = 0.35) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  labs(x = NULL, y = 'VAF') +
  ylim(0,1)
ggsave(paste(output_dir, '/figures/vaf_dotplot_all_smc56.eps', sep = ''),
       dpi = 320, width = 5, height = 3.7, units = 'in')