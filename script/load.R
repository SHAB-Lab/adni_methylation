library("limma")
library("minfi")
library("RColorBrewer")
library("missMethyl") # Can take a short time...
library("minfiData")
library("Gviz")
library("DMRcate")
library("DMRcatedata")
library("stringr")
library("mCSEA")

ann_epic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

data_dir <- "data/ADNI_iDAT_files"
list.files(data_dir, recursive = TRUE)

# READ DATA
targets <- read.metharray.sheet("data", pattern = "Sample_Sheet.csv")
# READ ONLY SMALL SET
tar_small <- targets[1:10, ]
rg_set <- read.metharray.exp(targets = tar_small)

# SUMMARIZE DATA
rg_set
pData(rg_set)
getManifest(rg_set)
