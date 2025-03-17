
# Cargar las librerias
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

# Lectura de datos
targets <- read.metharray.sheet("data", pattern = "Sample_Sheet.csv")

# Leer solo un pequeño set
tar_small <- targets[1:10, ]
rgSet <- read.metharray.exp(targets = tar_small)

# Resumen de datos
rgSet
pData(rgSet)
getManifest(rgSet)

# Para transformar los ensayos rojos y verdes en señales metiladas 
# y no metiladas respectivamente se utilizará el siguiente código:
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

# Resumir los valores de p en un punto único para simplificar
# la comparación entre muestras
# Examinar el promedio de los valores de p en todas las muestras
# para identificar cualquier muestra fallida
par(mar=c(9, 6, 4, 2))
barplot(colMeans(detP), las=2, cex.names=0.8)
mtext("Mean detection p-values", 2, 5)

# Para saber la razón exacta de por qué una muestra fallase el control de 
# calidad se puede usar la siguiente función: 
controlStripPlot(rgSet, controls="BISULFITE CONVERSION II") 

# Los gráficos de las diferentes sondas de control pueden ser exportadas como un docmuento pdf usando la función qcReport 
qcReport(rgSet, pdf= "qcReport.pdf") 


# Seleccionar las muestras que se mantendrán para un futuro análisis 
keep <- !colnames(rgSet) == "birth.11" 

# Crear el subconjunto de rgSet
rgSet <- rgSet[,keep] 

# Chequear que las muestras de mala calidad hayan sido removidas viendo el 
# número de columnas restantes
# Crear un subconjunto de los "target" también
targets <- targets[keep,] 

# Normalización de datos, los resultados quedarán en un objeto GenomicRatioSet 
mSetSq <- preprocessFunnorm(rgSet)

# Comparar los datos no normalizados para visualizar el efecto de la 
# normalización 
par(mfrow=c(1,2))

# Graficar distribuciones antes de la normalización para la muestra 1
plotBetasByType(MSet[,1],main="Raw")

# El objeto normalizado es un GenomicRatioSet que no contiene
# la información de la sonda necesaria, primero debemos extraerla del MethylSet.
typeI <- getProbeInfo(MSet, type = "I")[, c("Name","nCpG")]
typeII <- getProbeInfo(MSet, type = "II")[, c("Name","nCpG")]
probeTypes <- rbind(typeI, typeII)
probeTypes$Type <- rep(x = c("I", "II"), times = c(nrow(typeI), nrow(typeII)))

# Para graficar las distribuciones de los datos normalizados para la muestra 1
plotBetasByType(getBeta(mSetSq)[,1], probeTypes = probeTypes, main="Normalized",)

