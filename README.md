# Investigating potential effects of SMC5/6 genetic alterations in cancers with R
The following information is associated with the manuscript: 

[Disruption of the Structural Maintenance of Chromosomes 5/6 complex enables tumor mutagenesis](https://www.medrxiv.org/content/10.64898/2025.12.04.25341651v1.full)

Thi Tran<sup>1,2</sup>, Jiayi Fan<sup>3</sup>, Xiaolan Zhao<sup>3</sup>, Abby M. Green<sup>1,2</sup>

<sup>1</sup>Department of Pediatrics, Washington University School of Medicine, St. Louis, MO, USA

<sup>2</sup>Center for Genome Integrity, Siteman Cancer Center, Washington University School of Medicine, St. Louis, MO, USA

<sup>3</sup>Department of Molecular Biology, Memorial Sloan Kettering Cancer Center, New York, NY, USA

## System Requirement 
All the analysis have been run on a server with 12 CPUs and 128GB of RAM.

## Data 
* Somatic variant calls and copy number information of the TCGA projects were obtained using R package TCGA biolinks (version 2.32.0).
* Whole genome sequencing data of the International Cancer Genome Consortium (ICGC) was downloaded at <https://dcc.icgc.org/releases/PCAWG> (release 28).
* Clinical information, aneuploidy and ploidy scores for both TCGA and PCAWG (ICGC) dataset were obtained from cBioPortal on January 15th, 2024. 
* Somatic variant calling and clinical data of colorectal metastases were performed and reported by the Hartwig Medical Foundation (study number NCT01855477).
* To comply with data sharing agreement across different databases included in this study, we will not provide neither raw nor secondary datasets but will include R scripts used to generate results seen in the manuscript. 
