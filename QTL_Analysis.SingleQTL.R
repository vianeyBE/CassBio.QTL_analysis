# QTL Analysis using Rqtl package
# Author: Vianey Barrera (vpbarrerae@gmail.com / v.barrera@cgiar.org)
# The script runs the single QTL analysis (scanone) using data in mapqtl format.
# It computes QTLs, its Lod/Bayes intervals, scores and effects
# arguments:
# dircross: dir were data is located. It can be an empty string if the cross files include the path and name
# dir: dir to save results.
# dirfun: dir were functions for plot are located
# locfile: .loc file name, it can include the path
# mapfile: .map file name, it can include the path
# phenofile: .qua file name, it can include the path
# prefixResults: Results will have this prefix
# ncores: number of cores for permutation test
# step <- 0.5 : scanone arguments
# off.end <- 0 : scanone arguments
# error.prob <- 0.001 : scanone arguments
# alpha <- 0.1 : scanone arguments
# map.function <- "kosambi" : scanone arguments
# stepwidth <- "fixed" : scanone arguments
# model_scanone <- "normal" : scanone arguments
# n.perm <- 1000 : scanone arguments - For permutation test



single_qtl <- function(dircross, dir, dirfun, locfile, mapfile, phenofile, prefixResults, 
                       ncores, step=0.5, off.end=0, error.prob=0.001, alpha=0.1, n.perm=1000,
                       map.function="kosambi", stepwidth="fixed", model_scanone="normal"){
  
  # Loading Packages and variables ------------------------------------------
  
  library("qtl")
  library("tidyverse")
  library("snow")
  
  options(warn=-1)
  
  # Reading Cross  ----------------------------------------------------------
  
  cross <- read.cross(format = "mapqtl", 
                      dir = dircross, 
                      genfile = locfile,
                      mapfile = mapfile,
                      phefile = phenofile, 
                      genotypes = NULL)
  
  message("Data was loaded successfully!")
  
  cross <- jittermap(cross)
  cross <- drop.nullmarkers(cross)
  
  nind <- nind(cross)
  message(paste0("Number of individuals: ", nind))
  npheno <- nphe(cross)
  message(paste0("Number of traits: ", npheno-1))
  nmarker <- totmar(cross)
  message(paste0("Number of SNPs: ", nmarker))
  nchro <- nchr(cross)
  message(paste0("Number of chromosomes: ", nchro))
  
  #ngeno = length(levels(as.factor(cross$geno[[1]][["data"]])))
  
  # Genetic map data
  data.map <- pull.map(cross, as.table=TRUE) %>% rownames_to_column(var='mar')
  
  # Phenotype data
  data.pheno <- cross[["pheno"]]
  phenotypes <- colnames(data.pheno)
  
  
  # Interval Mapping --------------------------------------------------------
  
  cross <- calc.genoprob(cross, map.function = map.function, stepwidth = stepwidth,
                         step = step, error.prob = error.prob, off.end = off.end)
  
  out <- scanone(cross, method = "em", model = model_scanone,
                 pheno.col = c(2:dim(cross$pheno)[2]))
  
  message("Interval Mapping: done!")
  
  # Permutation Test --------------------------------------------------------
  
  message(paste0("Number of interations for permutation test: ", n.perm))
  
  out.perm <- scanone(cross, method = "em", n.perm = n.perm, n.cluster = ncores, 
                      pheno.col = c(2:dim(cross$pheno)[2]))
  
  message(paste0(capture.output(summary(out.perm, alpha = c(0.10, 0.05, 0.01))), collapse = "\n"))
  
  message("Permutation Test: done!")
  
  
  # Finding Peaks -----------------------------------------------------------
  
  # threshold
  lod.threshold <- melt(as.data.frame(print(summary(out.perm, alpha = alpha))), 
                        value.name = "lod", variable.name = "phenotype")
  
  message("Finding Peaks...")
  
  # Estimate picks 
  peaks <- (summary(object = out, 
                    perms = out.perm, format ="allpheno", 
                    threshold = as.numeric(lod.threshold[,2]))) %>% 
    melt(id.vars = c("chr", "pos"), 
         value.name = "lod", variable.name = "phenotype") %>%
    merge(lod.threshold, by = "phenotype", suffixes = c("", ".threshold") ) %>% 
    subset(lod > lod.threshold)
  
  peaks$pheno_num <- peaks$phenotype
  
  levels(peaks$pheno_num) <- c(1:length(levels(peaks$phenotype)))
  
  message(paste0("Peaks detected: ", dim(peaks)[1] ))
  
  # Estimate Intervals ------------------------------------------------------
  
  message("Estimating Intervals using LOD and Bayes method...")
  
  colnames_intervals <- c("phenotype", "pheno_num", "chr", "lod", 
                          "start.marker", "start.pos", 
                          "max.marker", "max.pos", 
                          "end.marker", "end.pos")
  
  # Estimate Intervals - LOD method
  lodint <- matrix(ncol=10, nrow = 0)
  colnames(lodint) <- colnames_intervals
  
  for (i in 1:dim(peaks)[1]){
    lodint.single <- lodint(out, chr = as.character(peaks$chr[i]), 
                            lodcolumn=as.numeric(peaks$pheno_num[i]), 
                            expandtomarkers = TRUE, drop = 1)
    
    lodint.temp <- data.frame(phenotype = as.character(peaks$phenotype[i]),
                              pheno_num = as.numeric(peaks$pheno_num[i]), 
                              chr = as.character(peaks$chr[i]), 
                              lod = peaks$lod[i], 
                              start.marker = row.names(lodint.single)[1], 
                              start.pos = lodint.single[1,2], 
                              max.marker = row.names(lodint.single)[2], 
                              max.pos = lodint.single[2,2], 
                              end.marker = row.names(lodint.single)[3], 
                              end.pos = lodint.single[3,2])
    
    lodint <- rbind(lodint, lodint.temp)
  }
  
  lodint <- distinct(lodint, pheno_num, chr, start.pos, max.pos, end.pos, .keep_all = TRUE) 
  
  # Find markers
  for (i in 1:dim(lodint)[1]){
    lodint[i,7] <- find.marker(cross, chr = lodint[i,3], pos = lodint[i,8])
  }
  
  # Estimate Intervals - BAYES method
  bayesint <- matrix(ncol=10, nrow = 0)
  
  colnames(bayesint) <- colnames_intervals
  
  for (i in 1:dim(peaks)[1]){
    bayesint.single = bayesint(out, chr = as.character(peaks$chr[i]), 
                               lodcolumn = as.numeric(peaks$pheno_num[i]), 
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
  
  bayesint <- distinct(bayesint, pheno_num, chr, start.pos, max.pos, end.pos, .keep_all = TRUE) 
  
  # Find markers
  for (i in 1:dim(lodint)[1]){
    bayesint[i,7] <-  find.marker(cross, chr = bayesint[i,3], pos = bayesint[i,8])
  }
  
  
  # QTL Effects -------------------------------------------------------------
  
  message("Computing QTL Effects")
  
  qtl.effects <-  cbind(lodint[,c(1,3,4,7)], additive = NA, dominance = NA, ratio = NA)
  
  for (i in 1:dim(lodint)[1]){
    
    mar <-  find.marker(cross, chr=lodint$chr[i], pos = lodint$max.pos[i])
    eff <-  effectplot(cross, mname1=mar, pheno.col = lodint$pheno_num[i]+1, draw = FALSE)
    qtl.effects$additive[i] <- (eff[["Means"]][1] - eff[["Means"]][4])/2
    qtl.effects$dominance[i] <- (eff[["Means"]][2] + eff[["Means"]][3])/2 - (eff[["Means"]][1] + eff[["Means"]][4])/2
    effectplot(cross, mname1=mar, pheno.col = lodint$pheno_num[i]+1)
    plotPXG(cross, marker = mar, pheno.col = lodint$pheno_num[i]+1)
    
  }
  
  qtl.effects$ratio = qtl.effects$dominance/abs(qtl.effects$additive)
  
  # Saving Data -------------------------------------------------------------

  outputdir <- paste0(dir, prefixResults,  "_Results/")
  outputplot <- paste0(dir, prefixResults,  "_Plots/")
  dir.create(outputdir, showWarnings = FALSE)
  dir.create(outputplot, showWarnings = FALSE)
  
  message("Saving Data...")
  
  # RData
  save.image(paste(prefixResults, format(Sys.Date(), "%m%d%Y"), "RData", sep = "."))
  
  
  # Save Files ####
  
  write.table(lodint,
              file = paste(outputdir, prefixResults, ".qtl.LodIntervals.txt", sep = ""), sep = "\t", 
              row.names = FALSE, quote = FALSE )
  write.table(bayesint,
              file = paste(outputdir, prefixResults, ".qtl.BayesIntervals.txt", sep = ""), sep = "\t", 
              row.names = FALSE, quote = FALSE )
  write.table(out,
              file = paste(outputdir, prefixResults, ".lod.scores.txt", sep = ""), sep = "\t", 
              row.names = TRUE, quote = FALSE )
  write.table(qtl.effects,
              file = paste(outputdir, prefixResults, ".qtl.effects.txt", sep = ""), sep = "\t", 
              row.names = FALSE, quote = FALSE )
    
  # Plots -------------------------------------------------------------------
  
  message(paste0("Saving Plots at: ", outputdir, " and ", outputplot))
  
  source(paste0(dirfun, "QTL_Analysis.singleQTL_plot.R"))
  genoPlot.pdf(cross)
  missGenoPlot.pdf(cross)
  plotPhenotypes.pdf(data.pheno)
  plotLOD.pdf(out)
  
  source(paste0(dirfun, "QTL_Analysis.singleQTL_map_plot.R") )
  plot.single.QTL(lodint, "Lod")
  plot.single.QTL(bayesint, "Bayes")
  
  message("Done!")
  

}

