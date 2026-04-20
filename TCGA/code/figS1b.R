## ENVIRONMENT ----
library('ggplot2')
library('tidyverse')
input_dir <- 'TCGA'
output_dir <- 'TCGA'
theme_general <- theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'bottom',
    legend.title = element_blank()
  )
groups_color <- c('#FAA0A0', '#A52A2A', '#80B19B')
names(groups_color) <- c('syn_mutation', 'del_mutation', 'nonsyn_nondel_mutation')
smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
          'NSMCE4A', 'NSMCE4B', 'EID3',
          'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')
## INPUT ----
cohort_sample <- 
  list.files(path = input_dir, full.names = TRUE, pattern = 'cohort_for_fpkm.csv$') %>%
  purrr::map(read.csv, header = TRUE) 
colnames(cohort_sample[[2]]) <- c('COHORT', 'alteration_type', 'Tumor_Sample_Barcode',
                                  'Hugo_Symbol', 'project')
cohort_sample <- bind_rows(cohort_sample)
smc56_gene_size <- data.frame(gene = smc56_gene, 
                              start_position = c(70258978,17663812,27224994,125091679,29264989,29264989,
                                                 121957091,104303739,104303739,
                                                 94618669,94618669,100912963,100912963),
                              end_position = c(70354873,17800242,27268772,125367120,29269822,29269822,
                                               121975217,104305205,104305205,
                                               94739436,94739436,100965134,100965134),
                              largest_transcript_size = c(1101,1091,266,247,304,304,
                                                          385,333,333,
                                                          1058,1058,1186,1186)) %>%
  mutate(gene_size = end_position-start_position+1)
master_snv <- read.csv(paste(input_dir, '/data/gdc_snv.csv', sep = ''),
                       header = TRUE)
## MUTATION FREQUENCY ----
df <- filter(master_snv, Hugo_Symbol %in% smc56_gene) %>%
  select(all_of(c('Tumor_Sample_Barcode','Hugo_Symbol','Chromosome','HGVSc'))) %>%
  unique() %>%
  mutate(Hugo_Symbol = with(., ifelse(Hugo_Symbol=='NDNL2',
                                      'NSMCE3', ifelse(Hugo_Symbol=='EID3', 
                                                       'NSMCE4B', ifelse(Hugo_Symbol=='ANKRD32',
                                                                         'SLF1', ifelse(Hugo_Symbol=='FAM178A',
                                                                                        'SLF2',Hugo_Symbol)))))) %>%
  inner_join(smc56_gene_size, by = join_by(Hugo_Symbol == gene), relationship = 'many-to-one') %>%
  select(all_of(c('Hugo_Symbol','Chromosome','HGVSc',
                  'start_position','end_position','gene_size','largest_transcript_size'))) %>%
  unique() %>%
  group_by(Hugo_Symbol,Chromosome,gene_size,largest_transcript_size) %>%
  summarise(n_snv = n()) %>%
  mutate(mutation_frequency=paste(format(round(100*n_snv/gene_size,3), nsmall = 3), '%', sep = '')) %>%
  select(all_of(c('Hugo_Symbol','Chromosome','largest_transcript_size',
                  'gene_size','n_snv','mutation_frequency')))

## GRAPH ----
df <- filter(cohort_sample, Hugo_Symbol %in% smc56_gene) %>%
  select(all_of(c('Tumor_Sample_Barcode','Hugo_Symbol','alteration_type'))) %>%
  unique() %>%
  summarise(n=n(),.by=c('Hugo_Symbol','alteration_type')) %>%
  mutate_at('Hugo_Symbol', ~ifelse(. == 'NDNL2', 'NSMCE3',
                                   ifelse(. == 'ANKRD32', 'SLF1',
                                          ifelse(. == 'FAM178A', 'SLF2',
                                                 ifelse(. == 'EID3', 'NSMCE4B', .))))) %>%
  mutate(Hugo_Symbol=fct_relevel(Hugo_Symbol,'SMC5','SMC6','NSMCE1','NSMCE2',
                                 'NSMCE3','NSMCE4A','NSMCE4B','SLF2','SLF1')) %>%
  mutate(alteration_type=fct_relevel(alteration_type,'syn_mutation','nonsyn_nondel_mutation','del_mutation'))
ggplot(df, aes(x=Hugo_Symbol,y=n,color=alteration_type)) +
  geom_bar(position='stack') +
  scale_fill_manual(values=groups_color) +
  theme_general +
  labs(y='Number of tumors',x=NULL)
ggsave(paste(output_dir, '/figures/separated_consequence.eps', sep = ''),
       dpi = 320)