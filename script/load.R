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
rgSet <- read.metharray.exp(targets = tar_small)

# SUMMARIZE DATA
rgSet
pData(rgSet)
getManifest(rgSet)

# Para transformar los ensayos rojos y verdes en señales metiladas 
# y no metiladas correspondientemente se utilizará el siguiente código:
MSet <- preprocessRaw(rgSet) 
MSet 
# Comparar con rgSet 

# Si corremos getMeth y getUnmeth podemos ver las intensidades de 
# las matrices metiladas y no metiladas.   

head(getMeth(MSet)[,1:3]) 
head(getUnmeth(MSet)[,1:3]) 

# Crear un RatioSet con ratioConvert 
ratioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
# Observar los cambios en el ensayo
ratioSet

gset <- mapToGenome(ratioSet)
gset

# Las funciones getBeta, getM y getCN permiten obtener los valores Beta, 
# M y el número de copias de la matrix 

beta <- getBeta(gset) 
head(beta) 
m <- getM(gset) 
head(m) 
cn <- getCN(gset) 
head(cn) 

#Control de Calidad 1
qc <- getQC(MSet) 
plotQC(qc) 

# Para calcular el valor de p
detP <- detectionP(rgSet) 
head(detP) 

# Resumir los valores de p en un puto único para simplificar
# la comparación entre muestras
# Examinar el promedio de los valores de p en todas las muestras
# para identificar cualquier muestra fallida
par(mar=c(9, 6, 4, 2))
barplot(colMeans(detP), las=2, cex.names=0.8)
mtext("Mean detection p-values", 2, 5)
