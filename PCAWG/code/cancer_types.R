# ENVIRONMENT ####
library('tidyverse')
library('ggplot2')
library('ggpubr',
        lib.loc = '/storage2/fs1/abby.green/Active/R_libraries/Bunny-Wunnies Freak Out.Bunny-Wunnies Freak Out')
working_dir <- '/storage2/fs1/abby.green/Active/Users/thi/SMC56_de_novo/2025/PCAWG'

# INPUT ####
# to keep the whitelist patients only
white_list_sample <- read.csv(paste(working_dir, '/data/white_list_sample_sheet.csv', sep = ''),
                              header = TRUE)
# to keep the primary samples only 
clinical_sample <- read.delim('/storage2/fs1/abby.green/Active/RawData/TCGA_PCAWG/cBioPortal/pancan_pcawg_2020/data_clinical_sample.txt',
                              header = TRUE, sep = '\t', comment.char = '#', fill = FALSE, quote = '')
primary_sample <- filter(clinical_sample, SAMPLE_TYPE=='Primary') %>%
  inner_join(white_list_sample, by = join_by(SAMPLE_ID==icgc_specimen_id)) %>%
  filter(is.na(PLOIDY) == FALSE)

smc56_smv <- read.csv(paste(working_dir, '/data/smc56_snv.csv', sep = ''),
                      header = TRUE) %>%
  mutate(COHORT = 'SmV') %>%
  filter(Tumor_Sample_Barcode %in% primary_sample$aliquot_id)

smc56_low_cnv <- read.csv(paste(working_dir, '/data/smc56_low_cnv.csv', sep = ''),
                          header = TRUE) %>%
  mutate(COHORT = 'low CNV') %>%
  filter(Tumor_Sample_Barcode %in% primary_sample$aliquot_id)

smc56_high_cnv <- read.csv(paste(working_dir, '/data/smc56_high_cnv.csv', sep = ''),
                           header = TRUE) %>%
  mutate(COHORT = 'high CNV') %>%
  filter(Tumor_Sample_Barcode %in% primary_sample$aliquot_id)
# ANALYSIS ####
smc56_wt <- filter(primary_sample, !aliquot_id %in% union(smc56_smv$Tumor_Sample_Barcode,
                                                          union(smc56_low_cnv$Tumor_Sample_Barcode,
                                                                smc56_high_cnv$Tumor_Sample_Barcode))) %>%
  mutate(COHORT = 'WT')
# no introns 
smv <- inner_join(primary_sample, smc56_smv,
                  by = join_by(aliquot_id == Tumor_Sample_Barcode)) %>%
  filter(Variant_Classification != 'Intron') %>%
  select(all_of(c('aliquot_id','CANCER_TYPE','COHORT'))) %>%
  unique()
low_cnv <- inner_join(primary_sample, smc56_low_cnv,
                    by = join_by(aliquot_id == Tumor_Sample_Barcode)) %>%
  select(all_of(c('aliquot_id','CANCER_TYPE','COHORT'))) %>%
  unique()
high_cnv <- inner_join(primary_sample, smc56_high_cnv,
                    by = join_by(aliquot_id == Tumor_Sample_Barcode)) %>%
  select(all_of(c('aliquot_id','CANCER_TYPE','COHORT'))) %>%
  unique()
wt <- select(smc56_wt, all_of(c('aliquot_id','CANCER_TYPE','COHORT'))) %>%
  unique()
df <- rbind(smv, wt) %>%
  rbind(low_cnv) %>%
  rbind(high_cnv) %>%
  mutate(facet = with(., ifelse(COHORT %in% c('low CNV', 'high CNV'),
                                'CNV', COHORT))) %>%
  group_by(CANCER_TYPE, COHORT, facet) %>%
  summarise(n = n()) %>%
  ungroup() %>% 
  mutate(CANCER_TYPE = forcats::fct_relevel(CANCER_TYPE, 'Hepatobiliary Cancer',
                                            'Pancreatic Cancer','Esophagogastric Cancer',
                                            'Breast Cancer','Prostate Cancer',
                                            'Mature B-cell lymphoma','Renal Cell Carcinoma',
                                            'Ovarian Cancer','Essential Thrombocythemia',
                                            'Intraductal Papillary Mucinous Neoplasm','Bone Cancer',
                                            'Embryonal Tumor','Head and Neck Cancer',
                                            'Mature B-Cell Neoplasms','Myelodysplastic/Myeloproliferative Neoplasms',
                                            'Non-Small Cell Lung Cancer','Colorectal Cancer',
                                            'Glioma', 'Lung Cancer',
                                            'Bladder Cancer','Soft Tissue Sarcoma',
                                            'Endometrial Cancer','Uterine Endometrioid Carcinoma',
                                            'Cervical Cancer','Medulloblastoma',
                                            'Thyroid Cancer',
                                            'Acute myeloid leukemia','Mastocytosis',
                                            'Posttransplant Lymphoproliferative Disorders',
                                            'Melanoma'))
color_group <- c('grey', '#A52A2A', '#C6B3DD', '#702963')
names(color_group) <- c('WT', 'SmV', 'low CNV', 'high CNV')
ggplot(df, aes(x=n, y=CANCER_TYPE, fill=COHORT)) +
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
  labs(x = 'Number of tumors', y = NULL) +
  facet_wrap(~factor(facet, levels = c('SmV','CNV', 'WT')), ncol = 3, scales = 'free_x') +
  scale_fill_manual(values = color_group)
ggsave(paste(working_dir, '/figures/no_intron_cancer_type.eps', sep = ''),
       dpi = 320, width = 6, height = 7, units = 'in')
