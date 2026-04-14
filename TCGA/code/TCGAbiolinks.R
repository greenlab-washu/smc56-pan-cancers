#This R script was to download TCGA data using TCGAbiolinks package 
#would be executed on the local machine but stored in storage platform 
#or run on storage platform if Docker image is available

if (!requireNamespace('BiocManager', quietly = TRUE))
	install.packages('BiocManager')
BiocManager::install('TCGAbiolinks')
#OR install from GitHub
install.packages('remotes')
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
install.packages('DT')

library('TCGAbiolinks')
library('tidyverse')
library('DT')

setwd('TCGA')

#only use data from TCGA project
tcga_projects <- getGDCprojects()$project_id %>% 
	stringr::str_extract('TCGA.+')
tcga_projects <- tcga_projects[!is.na(tcga_projects)] 
#expression data 
query_expression <- GDCquery(project = tcga_projects,
                             data.category = 'Transcriptome Profiling',
                             data.type = 'Gene Expression Quantification',
                             workflow.type = 'STAR - Counts',
                             access = 'open'
                             )
expression_id <- getResults(query_expression)
#single nucleotide variants 
query_snv <- GDCquery(project = tcga_projects, 
	                    data.category = 'Simple Nucleotide Variation',
	                    access = 'open'
                      )
snv_id <- getResults(query_snv)
#copy number variation 
query_cnv <- GDCquery(project = tcga_projects,
	                    data.category = 'Copy Number Variation',
	                    data.type = 'Gene Level Copy Number',
	                    access = 'open'
                      )
cnv_id <- getResults(query_cnv)

#only keep the sample that have all three of these info available 
common_patients <- intersect(
  intersect(cnv_id$cases.submitter_id, expression_id$cases.submitter_id),
  substr(snv_id$cases, 1, 12)
)
#need the query for common patients only before downloading data
query_expression <- GDCquery(project = tcga_projects,
                             data.category = 'Transcriptome Profiling',
                             data.type = 'Gene Expression Quantification',
                             access = 'open',
                             workflow.type = 'STAR - Counts',
                             barcode = common_patients)
expression_id <- getResults(query_expression) #need these dataframes to refer files' name 
                                              #and tumor/patient ID
write.csv(expression_id, 'expression_id.csv', row.names = FALSE)

query_snv <- GDCquery(project = tcga_projects,
	data.category = 'Simple Nucleotide Variation', 
	access = 'open',
	barcode = common_patients)
snv_id <- getResults(query_snv)
write.csv(snv_id, 'snv_id.csv', row.names = FALSE)

query_cnv <- GDCquery(project = tcga_projects, 
	data.category = 'Copy Number Variation', 
	data.type = 'Gene Level Copy Number',
	access = 'open',
	barcode = common_patients)
cnv_id <- getResults(query_cnv)
write.csv(cnv_id, 'cnv_id.csv', row.names = FALSE)

query_clinical <- GDCquery(project = tcga_projects,
                           data.category = 'Clinical',
                           access = 'open',
                           barcode = common_patients)
#download the data (SNV only - as of Jan 29th, 2024)
#download gene expression data
#create the manifest file needed for downloading data 
getManifest(query_experssion, save = TRUE) #need to go back and rename the file before proceeding
                                           #to prevent errors due to duplication 
getManifest(query_snv, save = TRUE)
getManifest(query_cnv, save = TRUE)
getManifest(query_clinical, save = TRUE)

# this line of code needs to be executed in the Terminal
gdc-client download -d TCGA/expression -m TCGA/gdc_manifest_expression.txt

# Clinical data for head-and-neck cancer 
query_clinical <- GDCquery(project = 'TCGA-HNSC',
                           data.category = 'Clinical')
getManifest(query_clinical, save=T)
