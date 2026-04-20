# ENVIRONMENT ####
library('tidyverse')
library('ggplot2')
working_dir <- 'TCGA'

# INPUT ####
snv_master <- read.csv(paste(working_dir, '/data/gdc_snv.csv', sep = ''), header = TRUE)
groups_color <- c('#A52A2A')
names(groups_color) <- c('del_mutation')

# ANALYSIS ####
driver_snv <- filter(snv_master, Hugo_Symbol %in% c('TP53', 'PIK3CA'))

driver_recurrent_snv <- driver_snv %>% 
  mutate(aa_position = stringr::str_extract(HGVSp_Short, '[:digit:]+')) %>%
  mutate(mutation_type = with(., ifelse(
    stringr::str_detect(HGVSp_Short, '=') == TRUE, 'syn_mutation', ifelse(
      IMPACT == 'HIGH', 'del_mutation', ifelse(
        SIFT == '', 'nonsyn_nondel_mutation', ifelse(
          stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') <= 0.05, 'del_mutation', 'nonsyn_nondel_mutation' 
        )
      )
    )))) %>%
  select(all_of(c('Hugo_Symbol', 'HGVSp_Short', 'aa_position', 'mutation_type', 'Tumor_Sample_Barcode'))) %>%
  group_by(Hugo_Symbol, HGVSp_Short, aa_position, mutation_type) %>%
  count()
write.csv(driver_recurrent_snv, paste(working_dir, '/data/driver_recurrent_mut.csv', sep = ''))

# PLOT ####
driver_recurrent_snv %>%
  filter(mutation_type=='del_mutation') %>%
  mutate_at('aa_position', as.numeric) %>%
  ggplot(aes(x = aa_position, y = n)) +
  geom_segment(aes(x = aa_position, xend = aa_position,
                   y = 0, yend = n, color = mutation_type)) +
  geom_point(aes(x = aa_position, y = n, color = mutation_type)) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_blank(),
    legend.position = 'bottom'
  ) +
  scale_color_manual(values = groups_color) +
  ylab('Number of occurrences') +
  xlab('Protein location') +
  facet_wrap(~factor(Hugo_Symbol,
                     levels = c('TP53', 'PIK3CA')), 
             scales = 'free', nrow = 2)
ggsave(paste(working_dir, '/figures/driver_genes_lollipop.eps', sep =''),
       dpi = 320)
