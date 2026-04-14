# August 28th, 2025 
# Thi Tran 
# identify patients with biallelic SMC5/6 mutations (complete loss of normal copy)
# focus on the PASS SMC5/6 variants
# what about intron mutations? 

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

smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
                'NSMCE4A', 'NSMCE4B', 'EID3',
                'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')
# INPUT ####
# do we have enough information to perform the analysis 
# 'observed VAF / tumor purity x local copy number = absolute number of chromatids 
# variant ploidy > local copy number - 0.5 is marked as biallelic' 
clinical <- read_tsv(paste(input_dir, '/metadata.tsv', sep = ''), na = c('null', '')) # this contains tumor purity 
post_biopsy <- read_tsv(paste(input_dir, '/post_biopsy_drugs.tsv', sep = ''), na = c('null', ''))
treatment <- read_tsv(paste(input_dir, '/treatment_responses.tsv', sep = ''), na = c('null', '')) 
# just to make sure that I dont take into account 
# patients without survival information 
condense_survival <- read.csv(paste(output_dir, '/data/survival.csv', sep = ''),
                              header = TRUE) %>% 
  filter(!is.na(time))
# CNV - region, not gene level 
cnv_files <- list.files(path = paste(input_dir, '/cnv', sep = ''),
                        pattern = 'purple.cnv.somatic.tsv', full.names = TRUE)
files_name <- cnv_files %>%
  gsub(paste(input_dir, '/cnv/', sep = ''), '', .) %>%
  gsub('\\.purple\\.cnv\\.somatic\\.tsv', '', .)
names(cnv_files) <- files_name
cnv_indv <- purrr::map(cnv_files, read.delim, header = TRUE, 
                       sep = '\t', comment.char = '#', fill = FALSE, quote = '')
cnv_master <- bind_rows(cnv_indv, .id = 'Tumor_Sample_Barcode') 
# SMC5/6 PASS SNV files
pass_snv_files <- list.files(paste(output_dir, '/data', sep = ''),
                        pattern = 'patients_pass_only', full.names = TRUE)
pass_snv_master<- purrr::map(pass_snv_files, read.csv, header = TRUE) %>%
  bind_rows(.id = NULL)
smc56_pass_snv <- filter(pass_snv_master, Hugo_Symbol %in% smc56_gene)
# TRANSFORM and ANALYSIS ----
smc56_mut_cn <- inner_join(smc56_pass_snv, cnv_master, 
                           by = join_by(Tumor_Sample_Barcode,Chromosome == chromosome),
                           relationship = 'many-to-many') %>%
  filter(start <= Start_Position) %>% # make sure that the SMC56 mut is within the CNV regions
  filter(Start_Position <= end)
# 'observed VAF / tumor purity x local copy number = absolute number of chromatids 
# variant ploidy > local copy number - 0.5 is marked as biallelic' 
df <- inner_join(smc56_mut_cn, clinical[,c('sampleId','tumorPurity')], by = join_by(Tumor_Sample_Barcode == sampleId)) %>%
  mutate(VAF = t_alt_count/t_depth) %>%
  mutate(variant_ploidy = (VAF/tumorPurity)*copyNumber) %>%
  mutate(biallelic = with(., ifelse(variant_ploidy > (copyNumber - .5),
                                    'yes','no'))) %>% 
  select(all_of(c('Tumor_Sample_Barcode','Hugo_Symbol','Chromosome','biallelic','HGVSp_Short',
                  'Variant_Classification','Variant_Type','Consequence',
                  'VARIANT_CLASS','Reference_Allele','Tumor_Seq_Allele2',
                  'VAF','SIFT','IMPACT'))) %>%
  unique()
filter(df, biallelic == 'yes') %>%
  write.csv(paste(output_dir, '/data/smc56_biallelic.csv', sep = ''), 
            row.names = FALSE)
