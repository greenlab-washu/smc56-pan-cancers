# highlighted regions and amino acid change annotations are Jiayi Fan's work
library('ggplot2')
library('tidyverse')
library('ggbreak')

input_dir <- 'TCGA'
output_dir <- 'TCGA'

## INPUT
groups_color <- c('#FAA0A0', '#A52A2A', '#80B19B')
names(groups_color) <- c('syn_mutation', 'del_mutation', 'nonsyn_nondel_mutation')
recurrent_mutations <- read.csv(paste(input_dir, '/separated_recurrent_smc56_mutations.csv', sep = '')) %>%
  select(-X) %>%
  mutate_at('Hugo_Symbol', ~ifelse(. == 'NDNL2', 'NSMCE3',
                                   ifelse(. == 'ANKRD32', 'SLF1',
                                          ifelse(. == 'FAM178A', 'SLF2',
                                                 ifelse(. == 'EID3', 'NSMCE4B', .)))))
recurrent_mutations %>%
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
  facet_wrap(~factor(Hugo_Symbol,
                     levels = c('SMC5', 'SMC6','NSMCE1', 'NSMCE2',
                                'NSMCE3', 'NSMCE4A', 'NSMCE4B',
                                'SLF1', 'SLF2')), scales = 'free_y', nrow = 9) 
ggsave(paste(output_dir, '/figures/smc56_lollipop.eps', sep = ''),
       dpi = 320, height = 11, width = 8.5, units = 'in')
#graph only SMC6 and break y-axis
recurrent_mutations %>% 
  filter(Hugo_Symbol == 'SMC6') %>%
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
  xlim(0, 1175) +
  scale_color_manual(values = groups_color) +
  ylab('number of occurrences') +
  xlab('Protein position') +
  scale_y_break(c(5, 25), scales = 0.15)
ggsave(paste(output_dir, '/figures/smc6_lollipop.eps', sep = ''),
       width = 8.5, height = 3, units = 'in', dpi = 320)
