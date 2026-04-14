# Oct 2nd, 2025 
# Build an oncoplot for PCAWG, just to satisfy my curiosity 
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

# read maf ----
df <- read.maf(maf = pass_smv_master, clinicalData = clinical)
plotmafSummary(df)
ggsave(paste(output_dir, '/figures/oncoplot_summary.eps', sep = ''), dpi = 320)