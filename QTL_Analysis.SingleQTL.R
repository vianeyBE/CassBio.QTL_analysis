library("qtl")

prefixResults <- "CM8996_2023" 
dir <- ""
#dir <- choose.dir(default = "", caption = "Select folder")
locfile <- gsub("\\\\", "/", file.choose()) 
mapfile <- gsub("\\\\", "/", file.choose()) 
phenofile <-gsub("\\\\", "/", file.choose()) 


cross = read.cross(format = "mapqtl", 
                 dir = dir, 
                 genfile = locfile,
                 mapfile = mapfile,
                 phefile = phenofile, 
                 genotypes = NULL)
cross = jittermap(cross)
cross = drop.nullmarkers(cross)

nind = nind(cross)
npheno = nphe(cross)
nmarker =totmar(cross)
nchro = nchr(cross)
ngeno = length(levels(as.factor(cross$geno[[1]][["data"]])))

# Genetic map data
data.map = pull.map(cross, as.table=TRUE)
data.map = cbind(mar = row.names(data.map), data.frame(data.map))
row.names(data.map) = c()

# Phenotype data
data.pheno = cross[["pheno"]]
phenotypes = colnames(data.pheno)

## Interval mapping - Standard Mapping ####
cross = calc.genoprob(cross, step = 0.5, error.prob = 0.001)
out = scanone(cross, method = "em", pheno.col=c(2:dim(cross$pheno)[2]))

## Permutation Test ####
out.perm <- scanone(cross, method = "em", n.perm = 1000,  pheno.col=c(2:dim(cross$pheno)[2]), n.cluster=ncores)
summary(out.perm, alpha=c(0.10, 0.05, 0.01)) 

# Stablish threshold 
alpha = 0.1
lod.threshold = melt(as.data.frame(print(summary(out.perm, alpha = alpha))), 
                     value.name = "lod", variable.name = "phenotype")

# Estimate picks 
peaks = (summary(object = out, 
                      perms = out.perm, 
                      format ="allpheno", threshold = as.numeric(lod.threshold[,2]))) %>% 
  melt(id.vars = c("chr", "pos"), value.name = "lod", variable.name = "phenotype") %>%
  merge(lod.threshold, by = "phenotype", suffixes = c("", ".threshold") ) %>% 
  subset(lod > lod.threshold)
peaks$pheno_num = peaks$phenotype
levels(peaks$pheno_num) = c(1:length(levels(peaks$phenotype)))

#Estimate Intervals - LOD method ####
lodint = matrix(ncol=10, nrow = 0)
colnames(lodint) =  c("phenotype", "pheno_num", "chr", "lod", 
                        "start.marker", "start.pos", 
                        "max.marker", "max.pos", 
                        "end.marker", "end.pos")

for (i in 1:dim(peaks)[1]){
  lodint.single = lodint(out, chr = as.character(peaks$chr[i]), lodcolumn=as.numeric(peaks$pheno_num[i]), 
                         expandtomarkers = TRUE, drop = 1)
  lodint.temp = data.frame(phenotype = as.character(peaks$phenotype[i]),
                      pheno_num = as.numeric(peaks$pheno_num[i]), 
                      chr = as.character(peaks$chr[i]), 
                      lod = peaks$lod[i], 
                      start.marker = row.names(lodint.single)[1], 
                      start.pos = lodint.single[1,2], 
                      max.marker = row.names(lodint.single)[2], 
                      max.pos = lodint.single[2,2], 
                      end.marker = row.names(lodint.single)[3], 
                      end.pos = lodint.single[3,2])
  lodint = rbind(lodint, lodint.temp)
}

lodint <- distinct(lodint,pheno_num,chr,start.pos,max.pos,end.pos,.keep_all = TRUE) 

# Find markers
for (i in 1:dim(lodint)[1]){
  lodint[i,7] = find.marker(cross, chr = lodint[i,3], pos = lodint[i,8])
}

#Estimate Intervals - BAYES method ####
bayesint = matrix(ncol=10, nrow = 0)
colnames(bayesint) =  c("phenotype", "pheno_num", "chr", "lod", 
                        "start.marker", "start.pos", 
                        "max.marker", "max.pos", 
                        "end.marker", "end.pos")

for (i in 1:dim(peaks)[1]){
  bayesint.single = bayesint(out, chr = as.character(peaks$chr[i]), lodcolumn=as.numeric(peaks$pheno_num[i]), 
                             expandtomarkers = TRUE)
  bayesint.temp = data.frame(phenotype = as.character(peaks$phenotype[i]),
                             pheno_num = as.numeric(peaks$pheno_num[i]), 
                             chr = as.character(peaks$chr[i]), 
                             lod = peaks$lod[i], 
                             start.marker = row.names(bayesint.single)[1], 
                             start.pos = bayesint.single[1,2], 
                             max.marker = row.names(bayesint.single)[2], 
                             max.pos = bayesint.single[2,2], 
                             end.marker = row.names(bayesint.single)[3], 
                             end.pos = bayesint.single[3,2])
  bayesint = rbind(bayesint, bayesint.temp)
}

bayesint <- distinct(bayesint,pheno_num,chr,start.pos,max.pos,end.pos,.keep_all = TRUE) 

# Find markers
for (i in 1:dim(lodint)[1]){
  bayesint[i,7] = find.marker(cross, chr = bayesint[i,3], pos = bayesint[i,8])
}

# QTL Effects ####
qtl.effects = cbind(lodint[,c(1,3,4,7)], additive = NA, dominance = NA, ratio = NA)

for (i in 1:dim(lodint)[1]){
  mar = find.marker(cross, chr=lodint$chr[i], pos = lodint$max.pos[i])
  eff = effectplot(cross, mname1=mar, pheno.col = lodint$pheno_num[i]+1, draw = FALSE)
  qtl.effects$additive[i] = (eff[["Means"]][1] - eff[["Means"]][4])/2
  qtl.effects$dominance[i] = (eff[["Means"]][2] + eff[["Means"]][3])/2 - (eff[["Means"]][1] + eff[["Means"]][4])/2
  effectplot(cross, mname1=mar, pheno.col = lodint$pheno_num[i]+1)
  plotPXG(cross, marker = mar, pheno.col = lodint$pheno_num[i]+1)
}

effectplot(cross, mname1="C07_8723674_C", pheno.col = 1+1)
plotPXG(cross, marker = "C07_13735570_G", pheno.col = 1+1)


qtl.effects$ratio = qtl.effects$dominance/abs(qtl.effects$additive)

#save.image(paste(prefixResults, ".RData", sep = ""))

# Save Files and Plots ####
dir.create(paste(dir, "results", sep = ""))
dir.create(paste(dir, "plots", sep = ""))

write.table(lodint, file = paste(dir, "results/", prefixResults, ".qtl.LodIntervals.txt", sep = ""), sep = "\t", 
            row.names = FALSE, quote = FALSE )
write.table(bayesint, file = paste(dir, "results/", prefixResults, ".qtl.BayesIntervals.txt", sep = ""), sep = "\t", 
            row.names = FALSE, quote = FALSE )
write.table(out, file = paste(dir, "results/", prefixResults, ".lod.scores.txt", sep = ""), sep = "\t", 
            row.names = TRUE, quote = FALSE )
write.table(qtl.effects, file = paste(dir, "results/", prefixResults, ".qtl.effects.txt", sep = ""), sep = "\t", 
            row.names = FALSE, quote = FALSE )


# Better to run in computer instead of server because some libraries are not possible to update in the server
source("qtl-plots-functions.R")
genoPlot.pdf(cross)
missGenoPlot.pdf(cross)
plotPhenotypes.pdf(data.pheno)
plotLOD.pdf(out)

source("qtl-map-plot.single-qtl.R") 
plot.single.QTL(lodint, "Lod")
plot.single.QTL(bayesint, "Bayes")


