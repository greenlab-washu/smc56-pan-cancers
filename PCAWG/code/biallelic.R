# Sep 8th, 2025 
# biallelic SMC5/6 SmV in PCAWG data 

# ENVIRONMENT ----
library('tidyverse')
library('ggplot2')
input_dir <- '/storage2/fs1/abby.green/Active/RawData/TCGA_PCAWG/PCAWG'
output_dir <- '/storage2/fs1/abby.green/Active/Users/thi/SMC56_de_novo/2025/PCAWG'
smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
                'NSMCE4A', 'NSMCE4B', 'EID3',
                'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')

# INPUT ----
white_list_sample <- read.csv(paste(output_dir, '/data/white_list_sample_sheet.csv', sep = ''),
                              header = TRUE)
## VAF ----
clinical_sample <- read.delim('/storage2/fs1/abby.green/Active/RawData/TCGA_PCAWG/cBioPortal/pancan_pcawg_2020/data_clinical_sample.txt', 
                              header = TRUE, sep = '\t', comment.char = '#', fill = FALSE, quote = '')
primary_sample <- filter(clinical_sample, SAMPLE_TYPE=='Primary') %>% 
  inner_join(white_list_sample, by = join_by(SAMPLE_ID==icgc_specimen_id)) %>% 
  filter(is.na(PLOIDY) == FALSE)

# this SmV master file contains all PASS SmV of primary tumors from white list patients 
snv_master <- read.csv(paste(output_dir, '/data/white_list_smv.csv', sep = ''), header = TRUE)
# the following just to documented how I came up with this SmV master file 
#raw_snv <- read.delim(paste(input_dir, '/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf', sep = ''), header = TRUE, sep = '\t', comment.char = '#', fill = FALSE, quote = '')
# only consider the primary samples from white list patients 
#snv_master <- filter(raw_snv, Tumor_Sample_Barcode %in% primary_sample$aliquot_id)
#write.csv(snv_master, paste(output_dir, '/data/white_list_smv.csv', sep = ''), row.names = FALSE)

## tumor purity ----
tumor_purity <- read.delim(paste(input_dir, '/consensus_cnv/consensus.20170217.purity.ploidy.txt', sep = ''),
                           header = TRUE, sep = '\t', comment.char = '#', fill = FALSE, quote = '')
## local copy number ----
# icgc 
icgc_files <- list.files(paste(input_dir, '/consensus_cnv/somatic_icgc', sep = ''),
                         full.names = TRUE)
filtered_icgc_files <- icgc_files[sapply(icgc_files, stringr::str_extract,
                                         '(?<=somatic_icgc/).+(?=\\.consensus\\.)') %in% primary_sample$aliquot_id]
# tcga
tcga_files <- list.files(paste(input_dir, '/consensus_cnv/somatic_tcga', sep = ''),
                         full.names = TRUE)
filtered_tcga_files <- tcga_files[sapply(tcga_files, stringr::str_extract,
                                         '(?<=somatic_tcga/).+(?=\\.consensus\\.)') %in% primary_sample$aliquot_id]
# combine them 
filtered_files <- combine(filtered_icgc_files, filtered_tcga_files)
files_name <- stringr::str_extract(
  filtered_files, '(?<=somatic_[:lower:]{4}/).+(?=\\.consensus\\.)'
)
names(filtered_files) <- files_name
rm(filtered_icgc_files, filtered_tcga_files, icgc_files, tcga_files)
# CNV - region wise, not gene wise in data_cleaning.R
cnv_indv <- purrr::map(filtered_files, read.delim, header = TRUE, 
                       sep = '\t', comment.char = '#', fill = FALSE, quote = '')
cnv_master <- bind_rows(cnv_indv, .id = 'Tumor_Sample_Barcode') 

# BIALLELIC selection ----
# 'observed VAF / tumor purity x local copy number = absolute number of chromatids 
# variant ploidy > local copy number - 0.5 is marked as biallelic'
# non-intronic mutations only 
smc56_smv <- filter(snv_master, Hugo_Symbol %in% smc56_gene)
smc56_mut_cn <- filter(smc56_smv, Variant_Classification != 'Intron') %>%
  inner_join(cnv_master, by = join_by(Tumor_Sample_Barcode,Chromosome == chromosome), 
             relationship = 'many-to-many') %>%
  filter(start <= Start_position) %>% # make sure that the SMC56 mut is within the CNV regions
  filter(Start_position <= end)
# 'observed VAF / tumor purity x local copy number = absolute number of chromatids 
# variant ploidy > local copy number - 0.5 is marked as biallelic' 
df <- inner_join(smc56_mut_cn, tumor_purity[,c('samplename','purity')], by = join_by(Tumor_Sample_Barcode == samplename)) %>%
  mutate(VAF = t_alt_count/(t_alt_count+t_ref_count)) %>%
  mutate(variant_ploidy = (VAF/purity)*total_cn) %>%
  mutate(biallelic = with(., ifelse(variant_ploidy > (total_cn - .5),
                                    'yes','no'))) %>% 
  select(all_of(c('Tumor_Sample_Barcode','Hugo_Symbol','Chromosome','biallelic',
                  'Variant_Classification','Variant_Type',
                  'Reference_Allele','Tumor_Seq_Allele2',
                  'Project_Code','Donor_ID'))) %>% # completely forget that we haven't accessed the predicted consequence of these mutations yet (VEP and SIFT score)
  unique()
filter(df, biallelic == 'yes') %>%
  write.csv(paste(output_dir, '/data/smc56_biallelic.csv', sep = ''), 
            row.names = FALSE)
# SURVIVAL ----