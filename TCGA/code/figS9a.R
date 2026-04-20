# ENVIRONMENT ----
library('MutationalPatterns')
library('BSgenome')
library('tidyverse')
library('NMF')
output_dir <- 'output'
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg38'
library(ref_genome, character.only = TRUE)
indel_cosmic_v3.2 <- get_known_signatures(source = 'COSMIC_v3.2', genome = 'GRCh37', muttype = 'indel')
# i know there's a mismatch here but there's only indels signatures for GRCh37
smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
                'NSMCE4A', 'NSMCE4B', 'EID3',
                'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')
color_code <- c('#D3E3CAFF','#BED6B3FF','#92A587FF',
                '#A1CAF6FF', '#6592D6FF',
                '#011C40FF',
                c(rep('#dddddc',12)))
names(color_code) <- c('ID7', #MMRd
                       'ID1','ID2', # replication slippage
                       'ID6','ID8', #defective non-homologous end-joining 
                       'ID3', #tobacco
                       'ID4','ID5','ID9','ID10',
                       'ID11','ID12','ID13','ID14',
                       'ID15','ID16','ID17','ID18')
# 3 groups ----
## INPUT ----
vcf_files <- list.files(paste(output_dir, '/TCGA/data', sep = ''),
                        pattern = 'cnv|wt|smv.vcf', full.names = T) 
cond <- gsub(paste(output_dir,'/TCGA/data/smc56_',sep=''),'',vcf_files) %>%
  gsub('.vcf','',.)
grl <- read_vcfs_as_granges(vcf_files,cond,ref_genome,
                            type = 'all',remove_duplicate_variants = F)
## ANALYSIS ----
indel_grl <- get_mut_type(grl, type = 'indel')
indel_context <- get_indel_context(indel_grl, ref_genome)
indel_counts <- count_indel_contexts(indel_context)
as.data.frame(indel_counts) %>%
  select(all_of(c('smv','cnv','wt'))) %>%
  plot_indel_contexts(condensed = TRUE)
ggsave(paste(output_dir,'/TCGA/figures/indel_count_snv_3groups.eps',sep=''),
       dpi=320, width = 10, height = 5, units = 'in')
fit <- fit_to_signatures(indel_counts, indel_cosmic_v3.2)

as.data.frame(fit$contribution) %>%
  mutate('id' = rownames(.)) %>%
  pivot_longer(!c(id), names_to = 'group', values_to = 'cont') %>%
  mutate(group = forcats::fct_relevel(group, 'smv','cnv','wt')) %>%
  mutate(id = forcats::fct_relevel(id, 'ID7', #MMRd
                                    'ID1','ID2', # replication slippage
                                    'ID6','ID8', #defective non-homologous end-joining 
                                    'ID3', #tobacco
                                    'ID4','ID5','ID9','ID10',
                                    'ID11','ID12','ID13','ID14',
                                    'ID15','ID16','ID17','ID18')) %>%
  ggplot(aes(x=group, y=cont, fill = id)) +
  geom_col(position = 'fill') +
  theme_bw() + 
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.title = element_blank(),
    text = element_text(family = 'Arial')
  ) +
  scale_fill_manual(values = color_code) +
  labs(x = NULL, y = 'Relative contribution', fill = 'IDsig')

ggsave(paste(output_dir, '/TCGA/figures/indel_sig_snv_3groups.eps',sep=''),
       dpi=320, width = 3.5, height=3.5, units='in')