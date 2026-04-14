# after downloading data using TCGAbiolinks.R
# based on un-normalized CNV 

library('tidyverse')


# INPUT ----

input_dir <- 'TCGA'
output_dir <- 'TCGA'

gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
          'NSMCE4A', 'NSMCE4B', 'EID3',
          'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')

#load a bunch of ID to identify samples 
cnv_id <- read.csv(paste(input_dir, '/cnv_id.csv', sep = '')) %>%
  filter(., stringr::str_detect(sample_type, 'Primary') == TRUE) 
#only keep the primary tumors
expression_id <- read.csv(paste(input_dir, '/expression_id.csv', sep = '')) %>%
  filter(., stringr::str_detect(sample_type, 'Primary') == TRUE)
#only keep the primary tumors
snv_id <- read.csv(paste(input_dir, '/snv_id.csv', sep = '')) %>%
  filter(., stringr::str_extract(cases, '(?<=TCGA.{9})[:digit:]{2}') %in% c('01', '03', '05', '09')) 
#only keep the primary tumors

#now come to common patients/tumors
common_patients <- intersect(
  intersect(cnv_id$cases.submitter_id, expression_id$cases.submitter_id),
  stringr::str_extract(snv_id$cases, '^.{12}')
)
common_tumors <- intersect(
  intersect(stringr::str_extract(cnv_id$sample.submitter_id, 'TCGA.{9}(01|03|05|09)[:upper:]'),
            expression_id$sample.submitter_id),
  stringr::str_extract(snv_id$cases, 'TCGA.{9}(01|03|05|09)[:upper:]')
)

# also need the aneuploidy information from cBioPortal
survival_files <- list.files(path = dir('cBioPortal/pan_can_atlas',full.names = TRUE),
                             full.names = TRUE, pattern = 'data_clinical_sample.txt')
survival_ploidy <- purrr::map(survival_files, read.delim, 
                              sep = '\t', header = TRUE, comment.char = '#') %>%
  bind_rows() %>% 
  filter(is.na(ANEUPLOIDY_SCORE) == FALSE)
common_patients <- intersect(common_patients, survival_ploidy$PATIENT_ID)

# SNV ----
snv_id <- snv_id %>%
  filter(., stringr::str_extract(cases, '^.{12}') %in% common_patients)

#load SNV input files
snv_files <- list.files(path = dir(paste(input_dir, '/snv', sep = ''),
                                   full.names = TRUE),
                        full.names = TRUE, pattern = '_masked.maf$')
files_name <- snv_files %>%
  gsub(paste(input_dir, '/snv/', sep = ''), '', .) %>%
  gsub('/.+', '', .)
names(snv_files) <- files_name
snv_files <- snv_files[names(snv_files) %in% snv_id$id]
snv_indv <- purrr::map(snv_files, read.delim, header = TRUE, 
                       sep = '\t', comment.char = '#', fill = FALSE, quote = '')

# filter samples with SMC56 mutations ----
chosen_column <- c(
  'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
  'Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2',
  'Match_Norm_Seq_Allele1','Match_Norm_Seq_Allele2',
  'Tumor_Validation_Allele1','Tumor_Validation_Allele2',
  'Match_Norm_Validation_Allele1','Match_Norm_Validation_Allele2',
  'Variant_Classification', 'Variant_Type', 'dbSNP_RS', 'Tumor_Sample_Barcode',
  'Mutation_Status', 'HGVSc', 'HGVSp', 'HGVSp_Short', 't_depth', 't_ref_count', 't_alt_count',
  'One_Consequence', 'Consequence', 'CLIN_SIG', 'BIOTYPE', 'CANONICAL', 'SIFT', 'PolyPhen',
  'PUBMED', 'IMPACT', 'VARIANT_CLASS', 'CONTEXT'
)
integer_column <- c(
  'Start_Position', 'End_Position', 't_depth', 't_ref_count', 't_alt_count'
)
charater_column <- c(
  'Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2',
  'Match_Norm_Seq_Allele1','Match_Norm_Seq_Allele2',
  'Tumor_Validation_Allele1','Tumor_Validation_Allele2',
  'Match_Norm_Validation_Allele1','Match_Norm_Validation_Allele2',
  'Hugo_Symbol', 'Chromosome', 'Variant_Classification', 'Variant_Type', 'dbSNP_RS',
  'Tumor_Sample_Barcode', 'Mutation_Status', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'One_Consequence',
  'Consequence', 'CLIN_SIG','BIOTYPE', 'CANONICAL', 'SIFT', 'PolyPhen', 'PUBMED', 'IMPACT', 
  'VARIANT_CLASS', 'CONTEXT'
)
snv_master <- 
  purrr::map(snv_indv, select, all_of(chosen_column)) %>%
  purrr::map(mutate_at, integer_column, as.integer) %>%
  purrr::map(mutate_at, charater_column, as.character) %>% 
  bind_rows(.id = NULL)
snv_master$Tumor_Sample_Barcode <- stringr::str_extract(
  snv_master$Tumor_Sample_Barcode, '^.{16}'
)
# select tumors based on their SMC5/6 mutations  
synSMC56 <- snv_master %>%
  filter(Hugo_Symbol %in% gene) %>% 
  filter(stringr::str_detect(HGVSp_Short, '=') == TRUE)
delSMC56 <- snv_master %>%
  filter(., Hugo_Symbol %in% gene) %>%
  filter(IMPACT == 'HIGH' |  (SIFT != '' &
    stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') <= 0.05))
nonsyn_nondelSMC <- snv_master %>%
  filter(Hugo_Symbol %in% gene) %>%
  filter(stringr::str_detect(HGVSp_Short, '=') == FALSE) %>%
  filter(IMPACT != 'HIGH' & 
           (stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') > 0.05 |
           SIFT == ''))
rm(charater_column, integer_column, chosen_column, snv_files)

# CNV ----
## load CNV input files ----
cnv_id <- cnv_id %>%
  filter(., cases.submitter_id %in% common_patients) %>%
  mutate('Tumor_Sample_Barcode' = 
           stringr::str_extract(sample.submitter_id, 'TCGA.{9}(01|03|05|09)[:upper:]'))

cnv_files <- list.files(path = dir(paste(input_dir, '/cnv', sep = ''),
                                   full.names = TRUE),
                        full.names = TRUE, pattern = 'gene_level_copy_number.v36.tsv')
files_name <- cnv_files %>%
  gsub(paste(input_dir, '/cnv/', sep = ''), '', .) %>%
  gsub('/.+', '', .)
names(cnv_files) <- files_name

new_cnv_id <- filter(cnv_id, id %in% names(cnv_files))
#since GDC only released ABSOLUTE files that were derived from
#alignment to hg19, ABSOLUTE files would be the least favored 
#(only use it if we have no other choice)
ASCAT3samples <- filter(new_cnv_id, analysis_workflow_type == 'ASCAT3')
ASCAT2samples <- filter(new_cnv_id, analysis_workflow_type == 'ASCAT2') %>%
  filter(., !Tumor_Sample_Barcode %in% ASCAT3samples$Tumor_Sample_Barcode)
ASCATsamples <- filter(new_cnv_id, analysis_workflow_type == 'AscatNSG') %>%
  filter(., !Tumor_Sample_Barcode %in% ASCAT2samples$Tumor_Sample_Barcode &
           !Tumor_Sample_Barcode %in% ASCAT3samples$Tumor_Sample_Barcode)
ABSOLUTEsamples <- filter(new_cnv_id, analysis_workflow_type == 'ABSOLUTE LiftOver') %>%
  filter(., !Tumor_Sample_Barcode %in% ASCATsamples$Tumor_Sample_Barcode &
           !Tumor_Sample_Barcode %in% ASCAT2samples$Tumor_Sample_Barcode &
           !Tumor_Sample_Barcode %in% ASCAT3samples$Tumor_Sample_Barcode)
cnv_files <- cnv_files[names(cnv_files) %in% c(ASCAT3samples$id, ASCAT2samples$id,
                                               ASCATsamples$id, ABSOLUTEsamples$id)]
new_cnv_id <- filter(cnv_id, id %in% names(cnv_files))
id <- new_cnv_id$id
names(id) <- new_cnv_id$Tumor_Sample_Barcode
names(cnv_files) <- names(id)
rm(ASCAT2samples, ASCAT3samples, ASCATsamples, ABSOLUTEsamples, id)
cnv_indv <- purrr::map(cnv_files, read.delim, header = TRUE, 
                       sep = '\t', comment.char = '#', fill = FALSE, quote = '')

cnv_master <- 
  purrr::map(cnv_indv, mutate_at, c('gene_id', 'gene_name', 'chromosome'), as.character) %>%
  purrr::map(mutate_at, c('start', 'end', 'copy_number', 'min_copy_number', 'max_copy_number'),
             as.integer) %>%
  bind_rows(.id = 'Tumor_Sample_Barcode')
rm(cnv_files, files_name)
## chose tumors based on their SMC5/6 copy number ----
lowSMC56 <- filter(cnv_master, gene_name %in% gene) %>%
  filter(is.na(copy_number) == TRUE | copy_number < 2)
highSMC56 <- filter(cnv_master, gene_name %in% gene) %>%
  filter(copy_number >= 3 & is.na(copy_number) == FALSE)

# OVERLAPPING ----
## syn - del -> del ----
# syn - nonsyn yet nondel -> nonsyn yet nondel 
syn_tumors <- unique(
  setdiff(
  synSMC56$Tumor_Sample_Barcode, delSMC56$Tumor_Sample_Barcode) %>%
  setdiff(nonsyn_nondelSMC$Tumor_Sample_Barcode))
## del - nonsyn yet nondel -> del ----
del_tumors <- unique(delSMC56$Tumor_Sample_Barcode)
del_high_tumors <- unique(
  intersect(
  delSMC56$Tumor_Sample_Barcode, highSMC56$Tumor_Sample_Barcode
))
## nonsyn yet nondel ----
nonsyn_nondel_tumors <- unique(
  setdiff(
  nonsyn_nondelSMC$Tumor_Sample_Barcode, delSMC56$Tumor_Sample_Barcode
))
## low CNA - high CNA -> omit for now (keep track) ----
low_high_tumors <- unique(
  intersect(
  lowSMC56$Tumor_Sample_Barcode, highSMC56$Tumor_Sample_Barcode
))
low_tumors <- unique(
  setdiff(
  lowSMC56$Tumor_Sample_Barcode, highSMC56$Tumor_Sample_Barcode
))
## high CNA -----
#high CNA - low CNA -> omit for now
high_tumors <- unique(
  setdiff(
  highSMC56$Tumor_Sample_Barcode, lowSMC56$Tumor_Sample_Barcode
))
## WT tumors - no alteration in SMC5/6 ----
wt_tumors <- unique(
  setdiff(
  snv_master$Tumor_Sample_Barcode, synSMC56$Tumor_Sample_Barcode
) %>%
  setdiff(delSMC56$Tumor_Sample_Barcode) %>%
  setdiff(nonsyn_nondelSMC$Tumor_Sample_Barcode) %>%
  setdiff(lowSMC56$Tumor_Sample_Barcode) %>%
  setdiff(highSMC56$Tumor_Sample_Barcode))

# READOUT ----

write.csv(snv_master, paste(output_dir, '/gdc_snv.csv', sep = ''),
          row.names = FALSE)

write.csv(cnv_master, paste(output_dir, '/gdc_cnv.csv', sep = ''),
          row.names = FALSE)

snv_id <-
  mutate(snv_id, 'Tumor_Sample_Barcode' = stringr::str_extract(cases, '^.{16}'))
## SYN ----
syn_cohort_tumors <- 
  filter(snv_master, Tumor_Sample_Barcode %in% syn_tumors) %>%
  mutate(COHORT = c(rep('syn', 52381))) %>%
  mutate(alteration_type = with(., ifelse(
    stringr::str_detect(HGVSp_Short, '=') == TRUE, 'syn_mutation', ifelse(
      IMPACT == 'HIGH', 'del_mutation', ifelse(
        SIFT == '', 'nonsyn_nondel_mutation', ifelse(
          stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') <= 0.05, 'del_mutation', 'nonsyn_nondel_mutation'
        )))
  ))) %>%
  select(all_of(c('COHORT', 'alteration_type', 'Tumor_Sample_Barcode', 'Hugo_Symbol'))) %>%
  inner_join(snv_id[, c('project', 'Tumor_Sample_Barcode')],
             by = join_by(Tumor_Sample_Barcode),
             relationship = 'many-to-many') %>%
  unique()
write.csv(syn_cohort_tumors, paste(output_dir, '/syn_cohort_for_fpkm.csv', sep = ''),
          row.names = FALSE)
## DEL ----
del_cohort_tumors <- 
  filter(snv_master, Tumor_Sample_Barcode %in% del_tumors) %>%
  mutate(COHORT = c(rep('del', 783374))) %>%
  mutate(alteration_type = with(., ifelse(
    stringr::str_detect(HGVSp_Short, '=') == TRUE, 'syn_mutation', ifelse(
      IMPACT == 'HIGH', 'del_mutation', ifelse(
        SIFT == '', 'nonsyn_nondel_mutation', ifelse(
          stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') <= 0.05, 'del_mutation', 'nonsyn_nondel_mutation'
        )))
  ))) %>%
  select(all_of(c('COHORT', 'alteration_type', 'Tumor_Sample_Barcode', 'Hugo_Symbol'))) %>%
  inner_join(snv_id[, c('project', 'Tumor_Sample_Barcode')],
             by = join_by(Tumor_Sample_Barcode),
             relationship = 'many-to-many') %>%
  unique()
write.csv(del_cohort_tumors, paste(output_dir, '/del_cohort_for_fpkm.csv', sep = ''),
          row.names = FALSE)
## NONSYN_NONDEL ----
nonsyn_nondel_cohort_tumors <- 
  filter(snv_master, Tumor_Sample_Barcode %in% nonsyn_nondel_tumors) %>%
  mutate(COHORT = c(rep('nonsyn_nondel', 113077))) %>%
  mutate(alteration_type = with(., ifelse(
    stringr::str_detect(HGVSp_Short, '=') == TRUE, 'syn_mutation', ifelse(
      IMPACT == 'HIGH', 'del_mutation', ifelse(
        SIFT == '', 'nonsyn_nondel_mutation', ifelse(
          stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') <= 0.05, 'del_mutation', 'nonsyn_nondel_mutation'
        )))
  ))) %>%
  select(all_of(c('COHORT', 'alteration_type', 'Tumor_Sample_Barcode', 'Hugo_Symbol'))) %>%
  inner_join(snv_id[, c('project', 'Tumor_Sample_Barcode')],
             by = join_by(Tumor_Sample_Barcode),
             relationship = 'many-to-many') %>%
  unique()
write.csv(nonsyn_nondel_cohort_tumors, paste(output_dir, '/nonsyn_nondel_cohort_for_fpkm.csv', sep = ''),
          row.names = FALSE)
## LOW_CNV -----
low_cnv_cohort_tumors <-
  filter(cnv_master, gene_name %in% gene) %>%
  filter(Tumor_Sample_Barcode %in% low_tumors) %>%
  mutate(COHORT = c(rep('low_cnv', 7812))) %>%
  mutate(alteration_type = with(., ifelse(
    is.na(copy_number) == TRUE, 'low_cna', ifelse(
      copy_number < 2, 'low_cna', ifelse(
        copy_number > 2, 'high_cna', 'normal_cna'
      ))
  ))) %>%
  select(all_of(c('COHORT', 'alteration_type', 'Tumor_Sample_Barcode', 'gene_name'))) %>%
  inner_join(cnv_id[, c('project', 'Tumor_Sample_Barcode')],
             by = join_by(Tumor_Sample_Barcode), relationship = 'many-to-many') %>%
  unique()
colnames(low_cnv_cohort_tumors) <- c('COHORT', 'alteration_type', 'Tumor_Sample_Barcode',
                                     'Hugo_Symbol', 'project')
write.csv(low_cnv_cohort_tumors, paste(output_dir, '/low_cnv_cohort_for_fpkm.csv',
                                       sep = ''), row.names = FALSE)
## HIGH_CNV ----
high_cnv_cohort_tumors <-
  filter(cnv_master, gene_name %in% gene) %>%
  filter(Tumor_Sample_Barcode %in% high_tumors) %>%
  mutate(COHORT = c(rep('high_cnv', 42408))) %>%
  mutate(alteration_type = with(., ifelse(
    is.na(copy_number) == TRUE, 'low_cna', ifelse(
      copy_number < 2, 'low_cna', ifelse(
        copy_number > 2, 'high_cna', 'normal_cna'
      ))
  ))) %>%
  select(all_of(c('COHORT', 'alteration_type', 'Tumor_Sample_Barcode', 'gene_name'))) %>%
  inner_join(cnv_id[, c('project', 'Tumor_Sample_Barcode')],
             by = join_by(Tumor_Sample_Barcode), relationship = 'many-to-many') %>%
  unique()
colnames(low_cnv_cohort_tumors) <- c('COHORT', 'alteration_type', 'Tumor_Sample_Barcode',
                                     'Hugo_Symbol', 'project')
write.csv(high_cnv_cohort_tumors, paste(output_dir, '/high_cnv_cohort_for_fpkm.csv',
                                       sep = ''), row.names = FALSE)

## WT -----
wt_cohort_tumors <- 
  filter(snv_master, Tumor_Sample_Barcode %in% wt_tumors) %>%
  mutate(COHORT = c(rep('wt', 262094))) %>%
  mutate(alteration_type = with(., ifelse(
    stringr::str_detect(HGVSp_Short, '=') == TRUE, 'syn_mutation', ifelse(
      IMPACT == 'HIGH', 'del_mutation', ifelse(
        SIFT == '', 'nonsyn_nondel_mutation', ifelse(
          stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') <= 0.05, 'del_mutation', 'nonsyn_nondel_mutation'
        )))
  ))) %>%
  select(all_of(c('COHORT', 'alteration_type', 'Tumor_Sample_Barcode', 'Hugo_Symbol'))) %>%
  inner_join(snv_id[, c('project', 'Tumor_Sample_Barcode')],
             by = join_by(Tumor_Sample_Barcode),
             relationship = 'many-to-many') %>%
  unique()
colnames(wt_cohort_tumors) <- c('COHORT', 'alteration_type', 'Tumor_Sample_Barcode',
                                     'Hugo_Symbol', 'project')
write.csv(wt_cohort_tumors, paste(output_dir, '/wt_cohort_for_fpkm.csv',
                                       sep = ''), row.names = FALSE)

## recurrent SMC5/6 mutations -----------------------------------------------
recurrent_mut <- filter(snv_master, Hugo_Symbol %in% gene) %>% 
  filter(stringr::str_extract(Tumor_Sample_Barcode, '.$') == 'A') %>%
  mutate(aa_position = stringr::str_extract(HGVSp_Short, '[:digit:]+')) %>%
  mutate(mutation_type = with(., ifelse(
    stringr::str_detect(HGVSp_Short, '=') == TRUE, 'syn_mutation', ifelse(
      IMPACT == 'HIGH', 'del_mutation', ifelse(
        SIFT == '', 'nonsyn_nondel_mutation', ifelse(
          stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') <= 0.05, 'del_mutation', 'nonsyn_nondel_mutation' 
        )
      )
    )))) %>%
  select(all_of(c('Hugo_Symbol', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'aa_position', 
                  'mutation_type', 'Tumor_Sample_Barcode'))) %>%
  group_by(Hugo_Symbol, HGVSc, HGVSp, HGVSp_Short, aa_position, mutation_type) %>%
  count()
write.csv(recurrent_mut, paste(output_dir, '/separated_recurrent_smc56_mutations.csv', sep = ''),
          row.names = FALSE)
