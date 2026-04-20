# Fig 3 - c, d, e, i
# Fig S8

library('tidyverse')
library('ggplot2')
library('ggpubr')
theme_general <- theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  )
input_dir <- 'TCGA'
output_dir <- 'TCGA'
smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
                'NSMCE4A', 'NSMCE4B', 'EID3',
                'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')
# INPUT ----
snv_master <- read.csv(paste(input_dir, '/data/gdc_snv.csv', sep = ''))
tmb_master <- read.csv(paste(input_dir,'/data/gdc_tmb.csv',sep=''), header = T)
cohort_sample <- 
  list.files(path = '/data',
             full.names = TRUE, pattern = 'cohort_for_fpkm.csv$') %>%
  purrr::map(read.csv, header = TRUE) 
colnames(cohort_sample[[2]]) <- colnames(cohort_sample[[4]])
cohort_master <- bind_rows(cohort_sample)

# distribution - all indels found ----
indel_all <- snv_master %>%
  mutate(status=with(.,ifelse(VARIANT_CLASS != 'SNV','contains indels','SNVs only'))) %>%
  select(all_of(c('Tumor_Sample_Barcode','status'))) %>%
  unique()
indel_smc56 <- filter(snv_master, 
                      Hugo_Symbol %in% smc56_gene & VARIANT_CLASS != 'SNV') %>%
  select(all_of(c('Tumor_Sample_Barcode'))) %>%
  mutate(status = 'contains SMC5/6 indels') %>%
  unique()
## 3i ----
df_plot<- tmb_master %>% 
  mutate(Tumor_Sample_Barcode = fct_reorder(Tumor_Sample_Barcode, -desc(TMB))) %>%
  mutate(quartile = with(., ifelse(TMB >= 4.833333e+00, 'last', 
                                   ifelse(TMB <= 1.000000e+00,
                                          'first','others')))) %>%
  left_join(unique(smc56_tumors[,c('Tumor_Sample_Barcode','smv_type')]), 
            by = join_by(Tumor_Sample_Barcode), relationship='one-to-many') %>%
  mutate(smv_type = replace_na(smv_type, 'No SMC5/6 SmV')) %>%
  mutate(status = with(.,ifelse(Tumor_Sample_Barcode %in% indel_smc56$Tumor_Sample_Barcode,
                                'SMC5/6 indel',str_replace(smv_type,'SmV','SNV')))) 
# STATS
unique(df_plot[,c('Tumor_Sample_Barcode','status','quartile')]) %>%
  summarise(.by=c('quartile','status'), n=n()) %>%
  filter(status != 'No SMC5/6 SNV')
unique(df_plot[,c('Tumor_Sample_Barcode','quartile')]) %>%
  summarise(.by=quartile, n=n())
# Chi-square test 
unique(df_plot[,c('Tumor_Sample_Barcode','TMB','quartile','status')]) %>%
  summarise(.by=c('quartile','status'), n=n()) %>%
  ungroup() %>%
  pivot_wider(names_from=status, values_from=n, values_fill=0) %>%
  select(-c('quartile')) %>%
  chisq.test()
# PLOT
df_plot %>%
  mutate(Tumor_Sample_Barcode = fct_reorder(Tumor_Sample_Barcode, -desc(TMB))) %>%
  mutate(status=fct_relevel(status,'Synonymous SMC5/6 SNV','Deleterious SMC5/6 SNV',
                            'NS/ND SMC5/6 SNV','SMC5/6 indel','No SMC5/6 SNV')) %>%
  ggplot(aes(x=Tumor_Sample_Barcode,y=TMB,color=status)) +
  geom_point() +
  geom_vline(xintercept = c('TCGA-XK-AAJR-01A','TCGA-09-2051-01A'),
             linetype = 'longdash') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = 'top') +
  scale_color_manual(values = c('#FAA0A0','#A52A2A','#80B19B','#333F48FF','#C8C9C7FF'))
ggsave(paste(output_dir, '/TCGA/figures/smc56indel_snv_tmb_all_tumors.eps',
             sep=''), dpi=320, width=7, height=3, units='in')

# distribution - SMC5/6 indels - number of SNV only ----
snv_count <- filter(snv_master, VARIANT_CLASS=='SNV') %>%
  select(all_of(c('Tumor_Sample_Barcode','Hugo_Symbol','Chromosome',
                  'Start_Position','End_Position','HGVSc'))) %>%
  unique() %>%
  summarise(.by=Tumor_Sample_Barcode, abs_snv_number=n()) %>%
  mutate(per_mb=abs_snv_number/30)
df_plot <- snv_count %>%
  mutate(status=with(.,ifelse(Tumor_Sample_Barcode %in% indel_smc56$Tumor_Sample_Barcode,
                              'contains SMC5/6 indels','others'))) %>%
  mutate(Tumor_Sample_Barcode = fct_reorder(Tumor_Sample_Barcode, -desc(per_mb))) %>%
  mutate(quartile = with(., ifelse(per_mb >= 4.633333e+00, 'last', 
                                   ifelse(per_mb <= 9.000000e-01,
                                          'first','others')))) 
## S8 ----
## STATS 
unique(df_plot[,c('Tumor_Sample_Barcode','status','quartile')]) %>%
  summarise(.by=c('quartile','status'),n=n())
unique(df_plot[,c('Tumor_Sample_Barcode','quartile')]) %>%
  summarise(.by=quartile, n=n())
# Chi-square test 
select(df_plot, all_of(c('Tumor_Sample_Barcode','per_mb','quartile','status'))) %>%
  unique() %>%
  summarise(.by=c('status','quartile'),n=n()) %>%
  ungroup() %>%
  pivot_wider(names_from=status, values_from=n, values_fill=0) %>%
  select(-c('quartile')) %>%
  chisq.test()
## PLOT 
filter(df_plot) %>%
  ggplot(aes(x=Tumor_Sample_Barcode,y=per_mb,color=status)) +
  geom_point() +
  geom_vline(xintercept = c('TCGA-B8-5552-01B','TCGA-AG-A00C-01A'),
             linetype = 'longdash') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = 'top') +
  scale_color_manual(values = c('#333E47','#C8C9C7FF'))
ggsave(paste(output_dir,'/TCGA/figures/smc56_indels_n_indels_all_tumors.eps',sep=''),
       dpi = 320, width = 7, height = 3, units = 'in')
# BOXPLOT ----
indel_count <- filter(snv_master, !VARIANT_CLASS %in% c('SNV','substitution') &
                        Reference_Allele != T & Tumor_Seq_Allele2 != T) %>% 
  select(all_of(c('Tumor_Sample_Barcode','Hugo_Symbol','Chromosome',
                  'Start_Position','End_Position','HGVSc'))) %>%
  unique() %>%
  summarise(.by=Tumor_Sample_Barcode, abs_indel_number=n()) %>%
  mutate(per_mb=abs_indel_number/30)
df_plot <- right_join(indel_count,unique(cohort_master[,c('COHORT','Tumor_Sample_Barcode')]),
                      by=join_by(Tumor_Sample_Barcode),relationship='one-to-many') %>%
  mutate_at('per_mb', replace_na, 0) 
## 3d - 3 groups ----
df_plot %>%
  mutate(status=with(.,ifelse(COHORT=='wt','Non-altered',
                              ifelse(COHORT %in% c('low_cnv','high_cnv'),'CNA','SmV')))) %>%
  mutate(status=fct_relevel(status,'SmV','CNA','Non-altered')) %>%
  select(all_of(c('status','Tumor_Sample_Barcode','per_mb'))) %>%
  unique() %>%
  ggplot(aes(x=status,y=log10(per_mb),color=status)) +
  geom_boxplot() +
  theme_general +
  stat_compare_means(comparisons = list(c('SmV','CNA'),c('CNA','Non-altered'),c('SmV','Non-altered')),
                     method = 'wilcox') +
  scale_color_manual(values=c('#A52A2A', '#702963', '#28446f'))
ggsave(paste(output_dir,'/TCGA/figures/log10indel_3groups.eps',sep=''),
       dpi=320, width=3.5, height=3.5, units='in')
## 3e - 4 groups ----
df_plot %>%
  filter(!COHORT %in% c('low_cnv','high_cnv')) %>%
  mutate(COHORT=fct_relevel(COHORT,'syn','del','nonsyn_nondel','wt')) %>%
  select(all_of(c('COHORT','Tumor_Sample_Barcode','per_mb'))) %>%
  unique() %>%
  ggplot(aes(x=COHORT,y=log10(per_mb),color=COHORT)) +
  geom_boxplot() +
  theme_general +
  stat_compare_means(comparisons = list(c('syn','del'),c('del','nonsyn_nondel'),
                                        c('del','wt'),
                                        c('syn','nonsyn_nondel'),c('nonsyn_nondel','wt'),
                                        c('syn','wt')),
                     method = 'wilcox') +
  scale_color_manual(values=c('#FAA0A0', '#A52A2A', '#80B19B','#28446f'))
ggsave(paste(output_dir,'/TCGA/figures/log10indel_4groups.eps',sep=''),
       dpi=320, width=4, height=3.5, units='in')
## 3c - log10TMB indels, syn, NS/ND, vs WT----
smc56_tumors <- filter(snv_master, Hugo_Symbol %in% smc56_gene) 

df <- inner_join(unique(cohort_master[,c('COHORT','Tumor_Sample_Barcode')]), tmb_master, 
                 by=join_by(Tumor_Sample_Barcode), relationship='many-to-one') %>%
  left_join(smc56_tumors, by=join_by(Tumor_Sample_Barcode), relationship='many-to-many') %>%
  mutate(status=with(.,ifelse(is.na(VARIANT_CLASS)==F,
                              ifelse(VARIANT_CLASS %in% c('deletion','insertion','indel'),
                                     'indel',COHORT),COHORT))) %>%
  select(all_of(c('Tumor_Sample_Barcode','status','TMB'))) %>%
  unique() 
df %>% 
  filter(status %in% c('del','indel','nonsyn_nondel','syn','wt')) %>%
  mutate(status=fct_relevel(status,'syn','del','nonsyn_nondel','indel','wt')) %>%
  ggplot(aes(x=status,y=log10(TMB),color=status)) +
  geom_boxplot() +
  theme_general +
  stat_compare_means(method = 'wilcox', 
                     comparisons = list(c('syn','indel'),c('del','indel'),
                                        c('nonsyn_nondel','indel'),c('indel','wt'))) +
  theme(legend.position = 'none') +
  scale_color_manual(values=c('#FAA0A0', '#A52A2A', '#80B19B','#333E47','#28446f'))
filter(df,status %in% c('syn','del','nonsyn_nondel','indel','wt')) %>%
  summarise(.by=status,n=n())
ggsave(paste(output_dir,'/TCGA/figures/log10tmb_indel_snv_wt.eps',sep=''),
        dpi = 320, width = 5, height = 4, units = 'in')
