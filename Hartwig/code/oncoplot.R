# Sep 16th, 2025 
# Build an oncoplot for Hartwig patients, including:
# types of treatments, specific treatments
# SMC5/6 alteration status 
# individual SBSsig 

# ENVIRONMENT ----
library('maftools')
library('tidyverse')
library('ggplot2')
input_dir <- '/storage2/fs1/abby.green/Active/Hartwig'
output_dir <- '/storage2/fs1/abby.green/Active/Users/thi/SMC56_de_novo/2025/Hartwig'
smc56_gene <- c('SMC5', 'SMC6', 'NSMCE1', 'NSMCE2', 'NSMCE3', 'NDNL2',
                'NSMCE4A', 'NSMCE4B', 'EID3',
                'SLF1', 'ANKRD32', 'SLF2', 'FAM178A')

# INPUT ----
# PASS only SmVs 
pass_smv_master <- list.files(paste(output_dir, '/data', sep = ''),
                              full.names = TRUE, pattern = 'patients_pass_only') %>%
  purrr::map(read.csv, header = TRUE) %>%
  bind_rows()
# clinical files 
treatment <- read_tsv(paste(input_dir, '/treatment_responses.tsv', sep = ''), na = c('null'))
clinical <- read_tsv(paste(input_dir, '/metadata.tsv', sep = ''), na = c('null')) %>%
  mutate(Tumor_Sample_Barcode = sampleId)
# cohort files 
smc56_del_master <- list.files(paste(output_dir, '/data', sep = ''),
                               pattern = 'patients.+del', full.names = TRUE) %>%
  map(read.csv, header = TRUE) %>%
  bind_rows() %>%
  mutate(smc56_COHORT = 'del')
smc56_nsnd_master <- list.files(paste(output_dir, '/data', sep = ''),
                                pattern = 'patients.+nsnd', full.names = TRUE) %>%
  map(read.csv, header = TRUE) %>%
  bind_rows() %>%
  mutate(smc56_COHORT = 'NS/ND')
smc56_syn <- read.csv(paste(output_dir, '/data/CPCT_patients_no_introns_syn_smc56.csv', sep = ''),
                      header = TRUE) %>%
  mutate(smc56_COHORT = 'syn')
smc56_smv_cohort <- bind_rows(unique(smc56_del_master[,c('Tumor_Sample_Barcode','smc56_COHORT')]),
                              unique(smc56_nsnd_master[,c('Tumor_Sample_Barcode','smc56_COHORT')])) %>%
  bind_rows(unique(smc56_syn[,c('Tumor_Sample_Barcode','smc56_COHORT')]))
# next i should try out this image() function 
# MAFTOOLS ----
df <- read.maf(maf = pass_smv_master, clinicalData = clinical)
plotmafSummary(df)
ggsave(paste(output_dir, '/figures/oncoplot_summary.eps', sep = ''), dpi = 320)

oncoplot(df, genes = smc56_gene)

image(x = 1:nrow(pws_mat_bg), y = 1:ncol(pws_mat_bg), z = pws_mat_bg, 
      axes = FALSE, xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
      col = "#ecf0f1")
image(x = 1:nrow(pws_mat), y = 1:ncol(pws_mat), z = pws_mat, 
      axes = FALSE, xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
      col = "#34495e", add = TRUE)
abline(h = (1:ncol(pws_mat)) + 0.5, col = "white")
abline(v = (1:nrow(pws_mat)) + 0.5, col = "white")
mtext(text = colnames(pws_mat), side = 2, at = 1:ncol(pws_mat), 
      line = 0.4, cex = fontSize, las = 2)
mtext(text = pw_pct, side = 4, at = 1:ncol(pws_mat), line = 0.4, 
      cex = fontSize, las = 2)
if (showTumorSampleBarcodes) {
  text(y = rep(0, nrow(pws_mat)), x = 1:nrow(pws_mat), 
       labels = rownames(pws_mat), srt = 45, font = 1, 
       cex = SampleNamefontSize, adj = 1, xpd = TRUE)