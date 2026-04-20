library('ggplot2')
library('tidyverse')

input_dir <- 'TCGA'
output_dir <- 'TCGA'

# INPUT ####
gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
          'NSMCE4A', 'NSMCE4B', 'EID3',
          'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')
group_color <- c('#FAA0A0', '#A52A2A', '#C6B3DD', '#702963')
names(group_color) <- c('syn', 'del', 'low_cnv', 'high_cnv')
#show the CNV and SNV cohorts separatedly
#load the input file
cohort_sample <- 
  list.files(path = input_dir, full.names = TRUE, pattern = 'cohort_for_fpkm.csv$') %>%
  purrr::map(read.csv, header = TRUE) 
colnames(cohort_sample[[2]]) <- colnames(cohort_sample[[4]])
cohort_sample <- bind_rows(cohort_sample)

# ANALYSIS ####
df1 <- filter(cohort_sample, COHORT %in% c('syn', 'del', 'low_cnv', 'high_cnv')) %>%
  mutate(facet = case_when(COHORT %in% c('syn', 'del') ~ 'SNA',
                           COHORT %in% c('low_cnv', 'high_cnv') ~ 'CNA')) 
df2 <- data.frame(project_id = c('TCGA-STAD', 'TCGA-GBM', 'TCGA-LUSC', 'TCGA-LUAD',
                                 'TCGA-BRCA', 'TCGA-BLCA', 'TCGA-UCEC', 'TCGA-COAD',
                                 'TCGA-SKCM', 'TCGA-LIHC', 'TCGA-READ', 'TCGA-HNSC',
                                 'TCGA-ACC', 'TCGA-OV', 'TCGA-ESCA', 'TCGA-KIRP',
                                 'TCGA-CESC', 'TCGA-PRAD', 'TCGA-LGG', 'TCGA-KIRC',
                                 'TCGA-DLBC', 'TCGA-SARC', 'TCGA-THYM', 'TCGA-PAAD',
                                 'TCGA-TGCT', 'TCGA-UCS', 'TCGA-PCPG', 'TCGA-UVM',
                                 'TCGA-KICH', 'TCGA-MESO', 'TCGA-LAML', 'TCGA-CHOL',
                                 'TCGA-THCA'),
                  primary_tissue = c('Stomach', 'Brain', 'Lung','Lung',
                                     'Breast', 'Bladder', 'Uterus', 'Colon', 
                                     'Skin', 'Liver', 'Rectum', 'Head and Neck', 
                                     'Adrenal gland', 'Ovary', 'Esophagus', 'Kidney', 
                                     'Cervix uteri', 'Prostate gland', 'Brain', 'Kidney', 
                                     'Lymphoma', 'Sacroma', 'Thymoma', 'Pancrea',
                                     'Testis', 'Uterus', 'Adrenal gland', 'Eye', 
                                     'Kidney', 'Mesothelium', 'AML', 'Bile ducts',
                                     'Thyroid gland'))
df <- inner_join(df1, df2, by = join_by(project == project_id), relationship = 'many-to-many') %>% 
  select(all_of(c('COHORT', 'primary_tissue', 'facet', 'Tumor_Sample_Barcode'))) %>%
  unique() %>%
  group_by(COHORT, primary_tissue, facet) %>%
  summarise(n = n()) 
#need to set 0 to some of the missing value in order to show the cancer in specified order
df0_syn <- data.frame(COHORT = rep('syn', 11),
                      primary_tissue = c('Lymphoma', 'Sacroma', 'Thymoma', 'Pancrea',
                                         'Testis', 'Adrenal gland', 'Eye',
                                         'Mesothelium', 'AML', 'Bile ducts', 'Thyroid gland'),
                      facet = rep('SNA', 11),
                      n = 0)
df0_del <- data.frame(COHORT = rep('del', 6),
                      primary_tissue = c('Lymphoma', 'Eye',
                                         'Mesothelium', 'AML', 'Bile ducts', 'Thyroid gland'),
                      facet = rep('SNA', 6),
                      n = 0)  
df <- bind_rows(df, df0_syn, df0_del)
df$COHORT <- factor(df$COHORT, levels = c('syn', 'del', 'low_cnv', 'high_cnv'))
df$primary_tissue <- factor(df$primary_tissue, levels = c('Uterus', 'Stomach', 'Lung',
                                                          'Colon', 'Bladder', 'Kidney', 'Breast',
                                                          'Cervix uteri', 'Esophagus', 'Liver', 'Rectum',
                                                          'Head and Neck', 'Prostate gland', 'Ovary',
                                                          'Brain', 'Skin', 'Adrenal gland', 'Sacroma', 'Pancrea', 
                                                          'Testis', 'Thymoma', 'Thyroid gland', 'Bile ducts',  
                                                          'AML', 'Eye', 'Lymphoma', 'Mesothelium'))
## GRAPH 
df %>% 
  ungroup() %>%
  ggplot(aes(x = primary_tissue, y = n, fill = COHORT)) +
  geom_col(position = position_dodge(preserve = 'single')) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        legend.position = 'bottom',
        legend.title = element_blank()
  ) +
  scale_fill_manual(values = group_color) +
  labs(x = NULL, y = 'Number of tumors') +
  facet_wrap(~factor(facet, levels = c('SNA', 'CNA')), nrow = 2, scales = 'free_y')
ggsave(paste(output_dir, '/figures/separated_abs_cnv_cancer_types.eps', sep = ''), dpi = 320)
