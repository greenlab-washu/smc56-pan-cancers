# August 20th, 2025 
# Thi Tran 
# SBSsig of patients with immunotherapy post treatment 

# ENVIRONMENT ####
library('MutationalPatterns')
library('BSgenome')
library('tidyverse')
library('ggplot2')
library('gtsummary')
library('ggpubr',
        lib.loc = '/storage2/fs1/abby.green/Active/R_libraries/Bunny-Wunnies Freak Out.Bunny-Wunnies Freak Out')
theme_general <- theme_classic() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

input_dir <- '/storage2/fs1/abby.green/Active/Hartwig'
output_dir <- '/storage2/fs1/abby.green/Active/Users/thi/SMC56_de_novo/2025/Hartwig'
#load reference genome - GRCh37, which is the reference genome used in Hartwig
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg19' 
library(ref_genome, character.only = TRUE)
# load COSMIC signature
snv_cosmic_3.2 <- get_known_signatures(muttype = 'snv', 'COSMIC_v3.2', genome = 'GRCh37')
indel_cosmic_3.2 <- get_known_signatures(muttype = 'indel', 'COSMIC_v3.2', genome = 'GRCh37')
## color ----
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
indel_color_code <- c('#362904FF',
                      '#E3B01DFF','#FEFB29FF',
                      '#BBD84EFF',
                      '#AFE8FCFF',
                      '#696C8FFF',
                      c(rep('#dddddc',12))
                      )
names(indel_color_code) <- c('ID7', # defective DNA mismatch repair
                             'ID1','ID2', # slippage during DNA replication of the replicated DNA strand
                             'ID6', # defective HR DNA damage repair
                             'ID8', # repair of DNA double strand breaks by NHEJ or mutations in TOP2A 
                             'ID17', # mutations in TOP2A
                             'ID3', # tobacco smoking
                             'ID13', # UV exposure
                             'ID18', # colibactin exposure
                             'ID4','ID5','ID9','ID10','ID11','ID12','ID14','ID15','ID16') # unknown
# INPUT ----
clinical <- read_tsv(paste(input_dir, '/metadata.tsv', sep = ''), na = c('null', ''))
post_biopsy <- read_tsv(paste(input_dir, '/post_biopsy_drugs.tsv', sep = ''), na = c('null', ''))
treatment <- read_tsv(paste(input_dir, '/treatment_responses.tsv', sep = ''), na = c('null', '')) 
# just to make sure that I dont take into account 
# patients without survival information 
condense_survival <- read.csv(paste(output_dir, '/data/survival.csv', sep = ''),
                              header = TRUE) %>% 
  filter(!is.na(time))
# cohort files
cohort_files <- list.files(paste(output_dir, '/data', sep = ''),
                           pattern = 'patients_.+_smc56.csv$', full.names = TRUE)
cohort_indv <- purrr::map(cohort_files, read.csv, header = TRUE)
# SmV
del_indv <- rbind(cohort_indv[[1]], cohort_indv[[5]]) %>%
  rbind(cohort_indv[[10]])
nsnd_indv <- rbind(cohort_indv[[2]], cohort_indv[[6]]) %>%
  rbind(cohort_indv[[11]])
syn_indv <- cohort_indv[[3]]
# CNV
low_cnv_indv <- cohort_indv[[9]]
high_cnv_indv <- cohort_indv[[8]]
# WT
wt_indv <- rbind(cohort_indv[[4]], cohort_indv[[7]]) %>%
  rbind(cohort_indv[[12]])
# load vcf files
vcfs <- list.files(paste(input_dir, '/snv', sep = ''), 
                        full.names = TRUE, pattern = 'purple.somatic.vcf$')
sample_names <- stringr::str_extract(vcfs, '(?<=snv/).+(?=\\.purple\\.somatic\\.vcf$)')
# load vcf as GRanges objects
grl <- read_vcfs_as_granges(vcfs, sample_names, ref_genome, type = 'all')
# ANALYSIS ----
## generate mutational matrix for each tumor from vcf ----
snv_grl <- get_mut_type(grl, type = 'snv')
mutation_matrix <- mut_matrix(vcf_list = snv_grl, ref_genome = ref_genome)
as.data.frame(mutation_matrix) %>%
  write.csv(paste(output_dir, '/data/indv_trinu_context.csv', sep = ''))
fit <- fit_to_signatures(mutation_matrix, snv_cosmic_3.2)
df <- apply(fit$contribution, 2, function(x) x / sum(x))
write.csv(df, paste(output_dir, '/data/sbssig_cont_indv_tumors.csv', sep = ''))
### no need for VCF anymore ----
df <- read.csv(paste(output_dir, '/data/sbssig_cont_indv_tumors.csv', sep = ''),
               header = TRUE)
# select those receiving immunotherapy post biopsy 
immuno_patients <- filter(post_biopsy, type == 'Immunotherapy') %>%
  filter(sampleId %in% condense_survival$sampleId)
df_immuno <- as.data.frame(df) %>% 
  select(any_of(c('X',immuno_patients$sampleId)))
df_plot <- as.data.frame(df_immuno) %>%
  pivot_longer(!c(X), names_to = 'sample', values_to = 'rel_cont') %>%
  mutate(group = with(., ifelse(sample %in% union(del_indv$Tumor_Sample_Barcode,
                                                  union(nsnd_indv$Tumor_Sample_Barcode, syn_indv$Tumor_Sample_Barcode)),
                                'SmV', 'others'))) %>%
  mutate(X = forcats::fct_relevel(X, 'SBS2','SBS13', #apobec
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

ggplot(df_plot, aes(x=sample, y=rel_cont, fill=X)) +
  geom_col(position='stack') + 
  theme_general +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = color_code) +
  labs(x='', y='Relative contribution') + 
  facet_wrap(~factor(group,levels = c('SmV', 'others')), scales = 'free_x')
ggsave(paste(output_dir, '/figures/sbssig_immuno_smv_others.eps', sep = ''),
       dpi = 320, width = 7, height = 5, units = 'in')

## heatmap ----
ggplot(df_plot,aes(x=X, y=sample, fill=rel_cont)) +
  geom_tile() +
  facet_grid(factor(group, levels = c('SmV', 'others'))~., 
             scales = 'free_y', space = 'free') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30),
        legend.position = 'bottom'
  ) +
  #guides(x=guide_axis(n.dodge = 2)) +
  scale_fill_distiller(name = 'Relative contribution', direction = -1, palette = 'GnBu')
ggsave(paste(output_dir, '/figures/sbssig_immuno_smv_others_heatmap.eps', sep = ''),
       dpi = 320, width = 8, height = 5.5, units = 'in')
## indel ----
indel_grl <- get_mut_type(grl, type = 'indel')
indel_context <- get_indel_context(indel_grl, ref_genome)
indel_counts <- count_indel_contexts(indel_context)
write.csv(indel_counts, paste(output_dir, '/data/indv_indel_count.csv', sep = ''))
fit <- fit_to_signatures(indel_counts, indel_cosmic_3.2)
df <- apply(fit$contribution, 2, function(x) x / sum(x))
write.csv(df, paste(output_dir, '/data/indel_sig_cont_indv_tumors.csv', sep = ''))
### no need for VCF anymore ----
indel_sig_master <- read.csv(paste(output_dir, '/data/indel_sig_cont_indv_tumors.csv', 
                                   sep = ''),
                             header = TRUE)
df_immuno <- as.data.frame(indel_sig_master) %>% 
  select(any_of(c('X',immuno_patients$sampleId)))
df_plot <- as.data.frame(df_immuno) %>%
  pivot_longer(!c(X), names_to = 'sample', values_to = 'rel_cont') %>%
  mutate(group = with(., ifelse(sample %in% union(del_indv$Tumor_Sample_Barcode,
                                                  union(nsnd_indv$Tumor_Sample_Barcode, syn_indv$Tumor_Sample_Barcode)),
                                'SmV', 'others'))) %>%
  mutate(X = forcats::fct_relevel(X,'ID7', # defective DNA mismatch repair
                                       'ID1','ID2', # slippage during DNA replication of the replicated DNA strand
                                       'ID6', # defective HR DNA damage repair
                                       'ID8', # repair of DNA double strand breaks by NHEJ or mutations in TOP2A 
                                       'ID17', # mutations in TOP2A
                                       'ID3', # tobacco smoking
                                       'ID13', # UV exposure
                                       'ID18', # colibactin exposure
                                       'ID4','ID5','ID9','ID10','ID11','ID12','ID14','ID15','ID16')) # unknown

ggplot(df_plot, aes(x=sample, y=rel_cont, fill=X)) +
  geom_col(position='stack') + 
  theme_general +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = indel_color_code) +
  labs(x='', y='Relative contribution') + 
  facet_wrap(~factor(group,levels = c('SmV', 'others')), scales = 'free_x')
ggsave(paste(output_dir, '/figures/indel_immuno_smv_others.eps', sep = ''),
       dpi = 320, width = 7, height = 5, units = 'in')
