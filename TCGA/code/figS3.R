library('tidyverse')
library('ggplot2')
library('ggbreak')

input_dir <- 'TCGA'
output_dir <- 'TCGA'

groups_color <- c('#FAA0A0', '#A52A2A', '#80B19B')
names(groups_color) <- c('syn_mutation', 'del_mutation', 'nonsyn_nondel_mutation')
cohesin_genes <- c('SMC1A', 'SMC1B', 'SMC3', 'STAG1', 'STAG2', 'STAG3',
                   'RAD21', 'NIPBL', 'MAU2', 'PDS5A', 'PDS5B', 'WAPL', 'CDCA5')
# INPUT ####
snv_master <- read.csv(paste(input_dir, '/gdc_snv.csv', sep = ''), header = TRUE)
# ANALYSIS ####
cohesin_snv <- filter(snv_master, Hugo_Symbol %in% cohesin_genes)
### Cohesin ----
cohesin_recurrent_snv <- cohesin_snv %>% 
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
write.csv(cohesin_recurrent_snv, paste(output_dir, '/data/cohesin_recurrent_mut.csv', sep = ''))

cohesin_recurrent_snv %>%
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
  ylab('number of occurrences') +
  xlab('Protein location') +
  facet_wrap(~factor(Hugo_Symbol,
                     levels = c('SMC1A', 'SMC1B', 'SMC3', 'STAG1', 'STAG2', 'STAG3',
                                'RAD21', 'NIPBL', 'MAU2', 'PDS5A', 'PDS5B', 'WAPL', 'CDCA5')), 
             scales = 'free', nrow = 7)
ggsave(paste(output_dir, '/figures/cohesin_lollipop.eps', sep = ''),
       dpi = 320, width = 11, height = 8.5, units = 'in')
#### plot NIPBL with y-axis break ----
cohesin_recurrent_snv %>% 
  filter(Hugo_Symbol == 'NIPBL') %>%
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
  ylab('number of occurrences') +
  xlab('Protein position') +
  scale_y_break(c(5, 25), scales = 0.15)
ggsave(paste(output_dir, '/figures/nipbl_lollipop.eps', sep = ''),
       dpi = 320, width = 5.5, height = 3, units = 'in')
