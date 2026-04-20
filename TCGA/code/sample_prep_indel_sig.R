# generate right input format to make S9a
# ENVIRONMENT ----
library('tidyverse')
vcf_header <- c(
  "##fileformat=VCFv4.2",
  "##source=myVCF",
  paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
)
output_dir <- 'output'
# INPUT ----
snv_master <- read.csv('TCGA/data/gdc_snv.csv',header = T)
cohort_sample <- 
  list.files(path = paste(output_dir, 'TCGA/data', sep = ''),
             full.names = TRUE, pattern = 'cohort_for_fpkm.csv$') %>%
  purrr::map(read.csv, header = TRUE) 
colnames(cohort_sample[[2]]) <- colnames(cohort_sample[[4]])
cohort_master <- bind_rows(cohort_sample)

# CONVERT ----
df <- filter(snv_master, Reference_Allele != T & Tumor_Seq_Allele2 != T) %>%
  mutate(POS=with(.,ifelse(VARIANT_CLASS %in% c('SNV','substitution'),Start_Position,Start_Position-1)),
         REF=with(.,ifelse(VARIANT_CLASS %in% c('SNV','substitution'),Reference_Allele,
                           ifelse(VARIANT_CLASS=='insertion',str_sub(CONTEXT,6,6),
                                  str_c(str_sub(CONTEXT,6,6),Reference_Allele)))),
         ALT=with(.,ifelse(VARIANT_CLASS %in% c('SNV','substitution'),Tumor_Seq_Allele2,
                           ifelse(VARIANT_CLASS=='deletion',str_sub(CONTEXT,6,6),
                                  str_c(str_sub(CONTEXT,6,6),Tumor_Seq_Allele2))))) %>%
  select(all_of(c('Tumor_Sample_Barcode','Chromosome','POS','REF','ALT'))) %>%
  mutate(ID='.', .after=POS) %>%
  mutate(
    Quality = ".",
    Filter = "PASS",
    Info = "."
  )
## 3 groups ----
smv_vcf <- filter(df,
                  Tumor_Sample_Barcode %in% union(cohort_sample[[1]]$Tumor_Sample_Barcode,
                                                  union(cohort_sample[[5]]$Tumor_Sample_Barcode,
                                                        cohort_sample[[4]]$Tumor_Sample_Barcode))) %>%
  select(-Tumor_Sample_Barcode)
writeLines(vcf_header, con = paste(output_dir, "/TCGA/data/smc56_smv.vcf", sep =''))
write.table(smv_vcf, paste(output_dir, "/TCGA/data/smc56_smv.vcf", sep = ''), 
            append = TRUE, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
cnv_vcf <- filter(df,
                  Tumor_Sample_Barcode %in% union(cohort_sample[[2]]$Tumor_Sample_Barcode,
                                                  cohort_sample[[3]]$Tumor_Sample_Barcode)) %>%
  select(-Tumor_Sample_Barcode)
writeLines(vcf_header, con = paste(output_dir, "/TCGA/data/smc56_cnv.vcf", sep =''))
write.table(cnv_vcf, paste(output_dir, "/TCGA/data/smc56_cnv.vcf", sep = ''), 
            append = TRUE, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
wt_vcf <- filter(df,
                  Tumor_Sample_Barcode %in% cohort_sample[[6]]$Tumor_Sample_Barcode) %>%
  select(-Tumor_Sample_Barcode)
writeLines(vcf_header, con = paste(output_dir, "/TCGA/data/smc56_wt.vcf", sep =''))
write.table(wt_vcf, paste(output_dir, "/TCGA/data/smc56_wt.vcf", sep = ''), 
            append = TRUE, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

## SMC5/6 SmVs vs WT ---- 
syn_vcf <- filter(df,
                  Tumor_Sample_Barcode %in% cohort_sample[[5]]$Tumor_Sample_Barcode) %>%
  select(-Tumor_Sample_Barcode)
writeLines(vcf_header, con = paste(output_dir, "/TCGA/data/smc56_syn.vcf", sep =''))
write.table(syn_vcf, paste(output_dir, "/TCGA/data/smc56_syn.vcf", sep = ''), 
            append = TRUE, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

del_vcf <- filter(df,
                  Tumor_Sample_Barcode %in% cohort_sample[[1]]$Tumor_Sample_Barcode) %>%
  select(-Tumor_Sample_Barcode)
writeLines(vcf_header, con = paste(output_dir, "/TCGA/data/smc56_del.vcf", sep =''))
write.table(del_vcf, paste(output_dir, "/TCGA/data/smc56_del.vcf", sep = ''), 
            append = TRUE, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

nsnd_vcf <- filter(df,
                  Tumor_Sample_Barcode %in% cohort_sample[[4]]$Tumor_Sample_Barcode) %>%
  select(-Tumor_Sample_Barcode)
writeLines(vcf_header, con = paste(output_dir, "/TCGA/data/smc56_nsnd.vcf", sep =''))
write.table(nsnd_vcf, paste(output_dir, "/TCGA/data/smc56_nsnd.vcf", sep = ''), 
            append = TRUE, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)