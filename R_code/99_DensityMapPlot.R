# 
# 
# 



# X: ---------------------------------------------------------------------------
dir <- 'D:/OneDrive - CGIAR/00_CassavaBioinformaticsPlatform/01_ACWP/03_F2_Phenotypic/06_LinkageMap/06_DensityMap'
data <- 'F2_AM1588_map.csv'
setwd(dir)

library(LinkageMapView)

mapData <- read.csv(data)

lmv.linkage.plot(mapData, 'F2_AM1588_DensityMap.pdf', denmap = T)
