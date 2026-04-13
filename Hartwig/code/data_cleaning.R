# ENVIRONMENT ####
library('tidyverse')
library('ggplot2')
input_dir <- '/storage2/fs1/abby.green/Active/Hartwig'
output_dir <- '/storage2/fs1/abby.green/Active/Users/thi/SMC56_de_novo/2025/Hartwig'
smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
          'NSMCE4A', 'NSMCE4B', 'EID3',
          'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')
# August 1st, 2025 - omit the patients with missing treatment response information ####
# load the survival files to only keep the patients having clinical info recorded
treatment <- read_tsv(paste(input_dir, '/treatment_responses.tsv', sep = ''), na = c('null')) %>%
  group_by(patientId) %>% 
  filter(any(!is.na(response))) %>%
  ungroup()
clinical <- read_tsv(paste(input_dir, '/metadata.tsv', sep = ''), na = c('null'))
## CNV ----
# load the CNV files 
cnv_files <- list.files(path = paste(input_dir, '/cnv', sep = ''),
                        full.names = TRUE)
filtered_cnv_files <- cnv_files[sapply(cnv_files, stringr::str_extract, 
                                       '(?<=cnv/).+(?=\\.purple\\.cnv\\.gene\\.tsv)') %in% treatment$sampleId]
files_name <- filtered_cnv_files %>%
  gsub(paste(input_dir, '/cnv/', sep = ''), '', .) %>%
  gsub('\\.purple\\.cnv\\.gene\\.tsv', '', .)
names(filtered_cnv_files) <- files_name
rm(cnv_files)
cnv_indv <- purrr::map(filtered_cnv_files, read.delim, header = TRUE, 
                       sep = '\t', comment.char = '#', fill = FALSE, quote = '')
cnv_master <- bind_rows(cnv_indv, .id = 'Tumor_Sample_Barcode') %>%
  select(all_of(c('Tumor_Sample_Barcode','chromosome','start','end','gene','minCopyNumber','maxCopyNumber',
                  'minRegionStartSupport','minRegionEndSupport','minMinorAlleleCopyNumber',
                  'minMinorAlleleCopyNumber'))) %>%
  unique()
smc56_cnv <- filter(cnv_master, gene %in% smc56_gene)
low_smc56_cnv <- filter(smc56_cnv, maxCopyNumber <= 1.5)
high_smc56_cnv <- filter(smc56_cnv, minCopyNumber >= 2.5)
## SmV ####
# load SNV input files
snv_files <- list.files(path = paste(input_dir, '/snv', sep = ''),
                        full.names = TRUE, pattern = 'maf$') 
filtered_snv_files <- snv_files[sapply(snv_files, stringr::str_extract, 
                                       '(?<=snv/).+(?=\\.snv\\.maf)') %in% treatment$sampleId]
files_name <- filtered_snv_files %>%
  gsub(paste(input_dir, '/snv/', sep = ''), '', .) %>%
  gsub('\\.snv\\.maf', '', .)
names(filtered_snv_files) <- files_name
rm(snv_files)

snv_indv1 <- purrr::map(filtered_snv_files[1:228], read.delim, header = TRUE, 
                       sep = '\t', comment.char = '#', fill = FALSE, quote = '')
snv_master1 <- snv_indv1 %>%
  purrr::map(mutate_at, 'PUBMED', as.character) %>%
  bind_rows %>%
  select(all_of(c('Hugo_Symbol','Chromosome','Start_Position','End_Position',
                  'Tumor_Sample_Barcode','HGVSp_Short',
                  'Variant_Classification','Variant_Type','Feature_type','Consequence','BIOTYPE','VARIANT_CLASS',
                  'Reference_Allele','Tumor_Seq_Allele2','t_depth','t_ref_count','t_alt_count',
                  'flanking_bps','HGVSc','FILTER','SIFT','IMPACT'))) %>%
  unique()
rm(snv_indv1)
gc()

snv_indv2 <- purrr::map(filtered_snv_files[229:400], read.delim, header = TRUE,
                        sep = '\t', comment.char = '#', fill = FALSE, quote = '')
snv_master2 <- snv_indv2 %>%
  purrr::map(mutate_at, 'PUBMED', as.character) %>%
  bind_rows %>%
  select(all_of(c('Hugo_Symbol','Chromosome','Start_Position','End_Position',
                  'Tumor_Sample_Barcode','HGVSp_Short',
                  'Variant_Classification','Variant_Type','Feature_type','Consequence','BIOTYPE','VARIANT_CLASS',
                  'Reference_Allele','Tumor_Seq_Allele2','t_depth','t_ref_count','t_alt_count',
                  'flanking_bps','HGVSc','FILTER','SIFT','IMPACT'))) %>%
  unique()
rm(snv_indv2)
gc()

snv_indv3 <- purrr::map(filtered_snv_files[400:456], read.delim, header = TRUE, 
                        sep = '\t', comment.char = '#', fill = FALSE, quote = '')
snv_master3 <- snv_indv3 %>%
  purrr::map(mutate_at, 'PUBMED', as.character) %>%
  bind_rows %>%
  select(all_of(c('Hugo_Symbol','Chromosome','Start_Position','End_Position',
                  'Tumor_Sample_Barcode','HGVSp_Short',
                  'Variant_Classification','Variant_Type','Feature_type','Consequence','BIOTYPE','VARIANT_CLASS',
                  'Reference_Allele','Tumor_Seq_Allele2','t_depth','t_ref_count','t_alt_count',
                  'flanking_bps','HGVSc','FILTER','SIFT','IMPACT'))) %>%
  unique()
rm(snv_indv3)
gc()

snv_master <- rbind(snv_master1, snv_master2) %>%
  rbind(snv_master3)
rm(snv_master1, snv_master2, snv_master3)
gc() # these steps hopefully help with too much memory required 
smc56_snv <- filter(snv_master, FILTER == 'PASS') %>%
  filter(Hugo_Symbol %in% smc56_gene)
# there are lots of SMC56 intron mutations, 
# we'll leave them out for now to increase the sample size of some cohorts
no_introns_smc56_snv <- filter(smc56_snv, Variant_Classification != 'Intron')
del_smc56_snv <- filter(no_introns_smc56_snv, 
                        IMPACT == 'HIGH' | stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') <= 0.05)
nsnd_smc56_snv <- filter(no_introns_smc56_snv,
                         !Tumor_Sample_Barcode %in% del_smc56_snv$Tumor_Sample_Barcode) %>%
  filter(stringr::str_detect(HGVSp_Short, '=')==FALSE) %>%
  filter(IMPACT != 'HIGH' &
           (stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') > 0.05 | 
           SIFT == ''))
syn_smc56_snv <- filter(no_introns_smc56_snv, stringr::str_detect(HGVSp_Short, '=')==TRUE &
                          !Tumor_Sample_Barcode %in% union(del_smc56_snv$Tumor_Sample_Barcode, 
                                                           nsnd_smc56_snv$Tumor_Sample_Barcode))
no_introns_smc56_wt <- filter(snv_master, !Tumor_Sample_Barcode %in% union(no_introns_smc56_snv$Tumor_Sample_Barcode,
                                                                union(low_smc56_cnv$Tumor_Sample_Barcode, high_smc56_cnv$Tumor_Sample_Barcode))) %>%
  mutate(COHORT = 'wt') %>%
  select(all_of(c('Tumor_Sample_Barcode', 'COHORT'))) %>%
  unique()
## OUTPUT ####
write.csv(cnv_master, paste(output_dir, '/data/cnv_master.csv', sep = ''),
          row.names = FALSE)
rm(cnv_master)
write.csv(low_smc56_cnv, paste(output_dir, '/data/low_cnv_smc56.csv', sep = ''),
          row.names = FALSE)
write.csv(high_smc56_cnv, paste(output_dir, '/data/high_cnv_smc56.csv', sep = ''),
          row.names = FALSE)
filter(snv_master, FILTER == 'PASS') %>%
  write.csv(paste(output_dir, '/data/pass_only_snv_master.csv', sep = ''),
            row.names = FALSE)
write.csv(snv_master, paste(output_dir, '/data/snv_master.csv', sep = ''),
          row.names = FALSE)
write.csv(nsnd_smc56_snv, paste(output_dir, '/data/no_introns_nsnd_smc56.csv', sep = ''),
          row.names = FALSE)
write.csv(del_smc56_snv, paste(output_dir, '/data/no_introns_del_smc56.csv', sep = ''),
          row.names = FALSE)
write.csv(syn_smc56_snv, paste(output_dir, '/data/no_introns_syn_smc56.csv', sep = ''),
          row.names = FALSE)
write.csv(smc56_snv, paste(output_dir, '/data/smc56_snv.csv', sep = ''),
          row.names = FALSE)
write.csv(no_introns_smc56_wt, paste(output_dir, '/data/no_introns_wt_smc56.csv', sep = ''),
          row.names = FALSE)
# August 13th - include all patients - deal with missing information later ####
treatment <- read_tsv(paste(input_dir, '/treatment_responses.tsv', sep = ''), na = c('null'))
clinical <- read_tsv(paste(input_dir, '/metadata.tsv', sep = ''), na = c('null'))
## CNV ----
# load the CNV files 
cnv_files <- list.files(path = paste(input_dir, '/cnv', sep = ''),
                        full.names = TRUE)
files_name <- cnv_files %>%
  gsub(paste(input_dir, '/cnv/', sep = ''), '', .) %>%
  gsub('\\.purple\\.cnv\\.gene\\.tsv', '', .)
names(cnv_files) <- files_name
cnv_indv <- purrr::map(cnv_files, read.delim, header = TRUE, 
                       sep = '\t', comment.char = '#', fill = FALSE, quote = '')
cnv_master <- bind_rows(cnv_indv, .id = 'Tumor_Sample_Barcode') %>%
  select(all_of(c('Tumor_Sample_Barcode','chromosome','start','end','gene','minCopyNumber','maxCopyNumber',
                  'minRegionStartSupport','minRegionEndSupport','minMinorAlleleCopyNumber',
                  'minMinorAlleleCopyNumber'))) %>%
  unique()
smc56_cnv <- filter(cnv_master, gene %in% smc56_gene)
low_smc56_cnv <- filter(smc56_cnv, maxCopyNumber <= 1.5)
high_smc56_cnv <- filter(smc56_cnv, minCopyNumber >= 2.5)
## SmV ####
# load SNV input files
snv_files <- list.files(path = paste(input_dir, '/snv', sep = ''),
                        full.names = TRUE, pattern = 'maf$') 
files_name <- snv_files %>%
  gsub(paste(input_dir, '/snv/', sep = ''), '', .) %>%
  gsub('\\.snv\\.maf', '', .)
names(snv_files) <- files_name
## ACTN and WIDE patients ----
snv_indv1 <- purrr::map(snv_files[1:26], read.delim, header = TRUE, # Patient ID starts with ACTN
                        sep = '\t', comment.char = '#', fill = FALSE, quote = '')
snv_master1 <- snv_indv1 %>%
  purrr::map(mutate_at, 'PUBMED', as.character) %>%
  bind_rows %>%
  select(all_of(c('Hugo_Symbol','Chromosome','Start_Position','End_Position',
                  'Tumor_Sample_Barcode','HGVSp_Short',
                  'Variant_Classification','Variant_Type','Feature_type','Consequence','BIOTYPE','VARIANT_CLASS',
                  'Reference_Allele','Tumor_Seq_Allele2','t_depth','t_ref_count','t_alt_count',
                  'flanking_bps','HGVSc','FILTER','SIFT','IMPACT'))) %>%
  unique()
rm(snv_indv1)
gc()

snv_indv4 <- purrr::map(snv_files[579:726], read.delim, header = TRUE, # patient ID starts with WIDE
                        sep = '\t', comment.char = '#', fill = FALSE, quote = '')
snv_master4 <- snv_indv4 %>%
  purrr::map(mutate_at, 'PUBMED', as.character) %>%
  bind_rows %>%
  select(all_of(c('Hugo_Symbol','Chromosome','Start_Position','End_Position',
                  'Tumor_Sample_Barcode','HGVSp_Short',
                  'Variant_Classification','Variant_Type','Feature_type','Consequence','BIOTYPE','VARIANT_CLASS',
                  'Reference_Allele','Tumor_Seq_Allele2','t_depth','t_ref_count','t_alt_count',
                  'flanking_bps','HGVSc','FILTER','SIFT','IMPACT'))) %>%
  unique()
rm(snv_indv4)
gc()

snv_master <- rbind(snv_master1, snv_master4) 
gc() # these steps hopefully help with too much memory required 
smc56_snv <- filter(snv_master, FILTER == 'PASS') %>%
  filter(Hugo_Symbol %in% smc56_gene)
# there are lots of SMC56 intron mutations, 
# we'll leave them out for now to increase the sample size of some cohorts
no_introns_smc56_snv <- filter(smc56_snv, Variant_Classification != 'Intron') #36 tumors
del_smc56_snv <- filter(no_introns_smc56_snv, 
                        IMPACT == 'HIGH' | stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') <= 0.05) #11 tumors
nsnd_smc56_snv <- filter(no_introns_smc56_snv,
                         !Tumor_Sample_Barcode %in% del_smc56_snv$Tumor_Sample_Barcode) %>%
  filter(stringr::str_detect(HGVSp_Short, '=')==FALSE) %>%
  filter(IMPACT != 'HIGH' &
           (stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') > 0.05 | 
              SIFT == '')) # 25 tumors 
syn_smc56_snv <- filter(no_introns_smc56_snv, stringr::str_detect(HGVSp_Short, '=')==TRUE &
                          !Tumor_Sample_Barcode %in% union(del_smc56_snv$Tumor_Sample_Barcode, 
                                                           nsnd_smc56_snv$Tumor_Sample_Barcode)) 
# hence there is 0 tumor left in the group of Syn 
no_introns_smc56_wt <- filter(snv_master, !Tumor_Sample_Barcode %in% union(no_introns_smc56_snv$Tumor_Sample_Barcode,
                                                                           union(low_smc56_cnv$Tumor_Sample_Barcode, high_smc56_cnv$Tumor_Sample_Barcode))) %>%
  mutate(COHORT = 'wt') %>%
  select(all_of(c('Tumor_Sample_Barcode', 'COHORT'))) %>%
  unique()

### OUTPUT ####
write.csv(cnv_master, paste(output_dir, '/data/all_patients_cnv_master.csv', sep = ''),
          row.names = FALSE)
write.csv(low_smc56_cnv, paste(output_dir, '/data/all_patients_low_cnv_smc56.csv', sep = ''),
          row.names = FALSE)
write.csv(high_smc56_cnv, paste(output_dir, '/data/all_patients_high_cnv_smc56.csv', sep = ''),
          row.names = FALSE)
filter(snv_master, FILTER == 'PASS') %>%
  write.csv(paste(output_dir, '/data/other_patients_pass_only_snv_master.csv', sep = ''),
            row.names = FALSE)
write.csv(snv_master, paste(output_dir, '/data/other_patients_snv_master.csv', sep = ''),
          row.names = FALSE)
write.csv(nsnd_smc56_snv, paste(output_dir, '/data/other_patients_no_introns_nsnd_smc56.csv', sep = ''),
          row.names = FALSE)
write.csv(del_smc56_snv, paste(output_dir, '/data/other_patients_no_introns_del_smc56.csv', sep = ''),
          row.names = FALSE)
#write.csv(syn_smc56_snv, paste(output_dir, '/data/other_patients_no_introns_syn_smc56.csv', sep = ''),
#          row.names = FALSE)
write.csv(smc56_snv, paste(output_dir, '/data/other_patients_smc56_snv.csv', sep = ''),
          row.names = FALSE)
write.csv(no_introns_smc56_wt, paste(output_dir, '/data/other_patients_no_introns_wt_smc56.csv', sep = ''),
          row.names = FALSE)

## CPCT patients ----
# because they require lots of computational resource and the system keeps crashing
snv_indv2 <- purrr::map(snv_files[27:474], read.delim, header = TRUE, # patient ID starts with CPCT
                        sep = '\t', comment.char = '#', fill = FALSE, quote = '')
snv_master2 <- snv_indv2 %>%
  purrr::map(mutate_at, 'PUBMED', as.character) %>%
  bind_rows %>%
  select(all_of(c('Hugo_Symbol','Chromosome','Start_Position','End_Position',
                  'Tumor_Sample_Barcode','HGVSp_Short',
                  'Variant_Classification','Variant_Type','Feature_type','Consequence','BIOTYPE','VARIANT_CLASS',
                  'Reference_Allele','Tumor_Seq_Allele2','t_depth','t_ref_count','t_alt_count',
                  'flanking_bps','HGVSc','FILTER','SIFT','IMPACT'))) %>%
  unique()
rm(snv_indv2)
gc()
smc56_snv <- filter(snv_master2, FILTER == 'PASS') %>%
  filter(Hugo_Symbol %in% smc56_gene)
# there are lots of SMC56 intron mutations, 
# we'll leave them out for now to increase the sample size of some cohorts
no_introns_smc56_snv <- filter(smc56_snv, Variant_Classification != 'Intron')
del_smc56_snv <- filter(no_introns_smc56_snv, 
                        IMPACT == 'HIGH' | stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') <= 0.05) # 17 tumors
nsnd_smc56_snv <- filter(no_introns_smc56_snv,
                         !Tumor_Sample_Barcode %in% del_smc56_snv$Tumor_Sample_Barcode) %>%
  filter(stringr::str_detect(HGVSp_Short, '=')==FALSE) %>%
  filter(IMPACT != 'HIGH' &
           (stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') > 0.05 | 
              SIFT == '')) # 103 tumors
syn_smc56_snv <- filter(no_introns_smc56_snv, stringr::str_detect(HGVSp_Short, '=')==TRUE &
                          !Tumor_Sample_Barcode %in% union(del_smc56_snv$Tumor_Sample_Barcode, 
                                                           nsnd_smc56_snv$Tumor_Sample_Barcode)) # 2 tumors
no_introns_smc56_wt <- filter(snv_master2, !Tumor_Sample_Barcode %in% union(no_introns_smc56_snv$Tumor_Sample_Barcode,
                                                                           union(low_smc56_cnv$Tumor_Sample_Barcode, high_smc56_cnv$Tumor_Sample_Barcode))) %>%
  mutate(COHORT = 'wt') %>%
  select(all_of(c('Tumor_Sample_Barcode', 'COHORT'))) %>%
  unique() # 13 tumors

### OUTPUT ----
filter(snv_master2, FILTER == 'PASS') %>%
  write.csv(paste(output_dir, '/data/CPCT_patients_pass_only_snv_master.csv', sep = ''),
            row.names = FALSE)
write.csv(snv_master2, paste(output_dir, '/data/CPCT_patients_snv_master.csv', sep = ''),
          row.names = FALSE)
write.csv(nsnd_smc56_snv, paste(output_dir, '/data/CPCT_patients_no_introns_nsnd_smc56.csv', sep = ''),
          row.names = FALSE)
write.csv(del_smc56_snv, paste(output_dir, '/data/CPCT_patients_no_introns_del_smc56.csv', sep = ''),
          row.names = FALSE)
write.csv(syn_smc56_snv, paste(output_dir, '/data/CPCT_patients_no_introns_syn_smc56.csv', sep = ''),
          row.names = FALSE)
write.csv(smc56_snv, paste(output_dir, '/data/CPCT_patients_smc56_snv.csv', sep = ''),
          row.names = FALSE)
write.csv(no_introns_smc56_wt, paste(output_dir, '/data/CPCT_patients_no_introns_wt_smc56.csv', sep = ''),
          row.names = FALSE)

## DRUP patients ----
# because they require lots of computational resource and the system keeps crashing
snv_indv3 <- purrr::map(snv_files[475:578], read.delim, header = TRUE, # patient ID starts with DRUP
                        sep = '\t', comment.char = '#', fill = FALSE, quote = '')
snv_master3 <- snv_indv3 %>%
  purrr::map(mutate_at, 'PUBMED', as.character) %>%
  bind_rows %>%
  select(all_of(c('Hugo_Symbol','Chromosome','Start_Position','End_Position',
                  'Tumor_Sample_Barcode','HGVSp_Short',
                  'Variant_Classification','Variant_Type','Feature_type','Consequence','BIOTYPE','VARIANT_CLASS',
                  'Reference_Allele','Tumor_Seq_Allele2','t_depth','t_ref_count','t_alt_count',
                  'flanking_bps','HGVSc','FILTER','SIFT','IMPACT'))) %>%
  unique()
rm(snv_indv3)
gc()
smc56_snv <- filter(snv_master3, FILTER == 'PASS') %>%
  filter(Hugo_Symbol %in% smc56_gene)
# there are lots of SMC56 intron mutations, 
# we'll leave them out for now to increase the sample size of some cohorts
no_introns_smc56_snv <- filter(smc56_snv, Variant_Classification != 'Intron')
del_smc56_snv <- filter(no_introns_smc56_snv, 
                        IMPACT == 'HIGH' | stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') <= 0.05) # 16 tumors
nsnd_smc56_snv <- filter(no_introns_smc56_snv,
                         !Tumor_Sample_Barcode %in% del_smc56_snv$Tumor_Sample_Barcode) %>%
  filter(stringr::str_detect(HGVSp_Short, '=')==FALSE) %>%
  filter(IMPACT != 'HIGH' &
           (stringr::str_extract(SIFT, '(?<=\\().+(?=\\))') > 0.05 | 
              SIFT == '')) # 49 tumors
# there are no tumors left in this SYN group
syn_smc56_snv <- filter(no_introns_smc56_snv, stringr::str_detect(HGVSp_Short, '=')==TRUE &
                          !Tumor_Sample_Barcode %in% union(del_smc56_snv$Tumor_Sample_Barcode, 
                                                           nsnd_smc56_snv$Tumor_Sample_Barcode))
no_introns_smc56_wt <- filter(snv_master3, !Tumor_Sample_Barcode %in% union(no_introns_smc56_snv$Tumor_Sample_Barcode,
                                                                            union(low_smc56_cnv$Tumor_Sample_Barcode, high_smc56_cnv$Tumor_Sample_Barcode))) %>%
  mutate(COHORT = 'wt') %>%
  select(all_of(c('Tumor_Sample_Barcode', 'COHORT'))) %>%
  unique() # 3 tumors

### OUTPUT ----
filter(snv_master3, FILTER == 'PASS') %>%
  write.csv(paste(output_dir, '/data/DRUP_patients_pass_only_snv_master3.csv', sep = ''),
            row.names = FALSE)
write.csv(snv_master3, paste(output_dir, '/data/DRUP_patients_snv_master3.csv', sep = ''),
          row.names = FALSE)
write.csv(nsnd_smc56_snv, paste(output_dir, '/data/DRUP_patients_no_introns_nsnd_smc56.csv', sep = ''),
          row.names = FALSE)
write.csv(del_smc56_snv, paste(output_dir, '/data/DRUP_patients_no_introns_del_smc56.csv', sep = ''),
          row.names = FALSE)
# there is no tumors with only SYN SMC5/6 SmV
#write.csv(syn_smc56_snv, paste(output_dir, '/data/DRUP_patients_no_introns_syn_smc56.csv', sep = ''),
#          row.names = FALSE)
write.csv(smc56_snv, paste(output_dir, '/data/DRUP_patients_smc56_snv.csv', sep = ''),
          row.names = FALSE)
write.csv(no_introns_smc56_wt, paste(output_dir, '/data/DRUP_patients_no_introns_wt_smc56.csv', sep = ''),
          row.names = FALSE)
