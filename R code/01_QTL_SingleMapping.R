# Short name: QTL Analysis using Rqtl package
# Description: It computes single QTLs, Lod/Bayes intervals, scores, and effects using a mapqtl format
# Output:
#
# Author: Vianey Barrera (vpbarrerae@gmail.com / v.barrera@cgiar.org)
# 
# Arguments:
# 1. Single QTL function:
# dir: Directory (dir) to save results
# dircross: dir where data is located. It can be empty if the cross files include the path and name.
# dirfun: dir were functions for plotting are located.
# locfile: .loc file name. It can include the path.
# mapfile: .map file name. It can include the path.
# phenofile: .qua file name. It can include the path.
# prefixResults: The prefix that the results will have.
# ncores: Number of cores for permutation test.
#
# 2 scanone function:
# step: Step size in cM used for interval mapping (default = 0.5).
# off.end: A value added to each end of the chromosome, to extend the genetic map (default = 0).
# error.prob: The probability of a genotyping error (default = 0.001).
# alpha: The significance level to use for peak detection (default = 0.1).
# n.perm: Number of permutations to use for significance testing (default = 1000).
# map.function: The genetic map function to use (default = "kosambi").
# stepwidth: Method used to compute step size. Can be "fixed" or "cov" (default = "fixed").
# model_scanone: Model used for single-QTL scan (Options: "normal", "binary" or "poisson". Default = "normal").



##### To do #####
# 1: Add heritability calculations



# 0: Function init -------------------------------------------------------------

single_qtl <- function(dir, dircross, dirfun, locfile, mapfile, phenofile, prefixResults, ncores,
                       step = 0.5, off.end = 0, error.prob = 0.001, alpha = 0.1, n.perm = 1000,
                       map.function = "kosambi", stepwidth = "fixed", model_scanone = "normal"){
  
  
  
  # 1: Loading packages and variables ------------------------------------------
  if (!require(qtl)) install.packages("qtl")
  if (!require(reshape2)) install.packages("reshape2")
  if (!require(tidyverse)) install.packages("tidyverse")
  if (!require(snow)) install.packages("snow")
  if (!require(viridis)) install.packages("viridis")
  
  library(qtl)
  library(reshape2)
  library(tidyverse)
  library(snow)
  library(viridis)
  
  options(warn = -1)
  
  
  
  # 2: Reading cross -----------------------------------------------------------
  
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
  npheno <- nphe(cross)
  nmarker <- totmar(cross)
  nchro <- nchr(cross)
  
  message(paste0("Number of individuals: ", nind, "\n",
                 "Number of traits: ", npheno - 1, "\n",
                 "Number of SNPs: ", nmarker, "\n",
                 "Number of chromosomes: ", nchro, "\n"))
  
  ngeno <- length(levels(as.factor(cross$geno[[1]][["data"]])))
  
  # Genetic map data
  data.map <- pull.map(cross, as.table = T) %>% rownames_to_column(var = 'mar')
  
  # Phenotype data
  data.pheno <- cross[["pheno"]]
  phenotypes <- colnames(data.pheno)
  
  
  
  # 3: Interval mapping --------------------------------------------------------
  
  cross <- calc.genoprob(cross, map.function = map.function, stepwidth = stepwidth,
                         step = step, error.prob = error.prob, off.end = off.end)
  
  out <- scanone(cross, method = "em", model = model_scanone, pheno.col = c(2:dim(cross$pheno)[2]))
  
  message("Interval mapping: Done!")
  
  
  
  # 4: Permutation test --------------------------------------------------------
  
  message(paste0("Number of interations for permutation test: ", n.perm))
  
  out.perm <- scanone(cross, method = "em", n.perm = n.perm, n.cluster = ncores, 
                      pheno.col = c(2:dim(cross$pheno)[2]))
  
  message(paste0(capture.output(summary(out.perm, alpha = c(0.10, 0.05, 0.01))), collapse = "\n"))
  
  message("Permutation test: Done!")
  
  
  
  # 5: Finding peaks -----------------------------------------------------------
  # Threshold
  lod.threshold <- melt(as.data.frame(print(summary(out.perm, alpha = alpha))), 
                        value.name = "lod", variable.name = "phenotype")
  
  message("Finding peaks...")
  
  # Estimate picks 
  peaks <- summary(object = out, perms = out.perm, format = "allpheno", 
                    threshold = as.numeric(lod.threshold[,2])) %>% 
    melt(id.vars = c("chr", "pos"), value.name = "lod", variable.name = "phenotype") %>%
    merge(lod.threshold, by = "phenotype", suffixes = c("", ".threshold")) %>% 
    subset(lod > lod.threshold)
  
  peaks$pheno_num <- peaks$phenotype
  
  levels(peaks$pheno_num) <- c(1:length(levels(peaks$phenotype)))
  
  message(paste0("Peaks detected: ", dim(peaks)[1] ))
  
  
  
  # 6: Estimate intervals ------------------------------------------------------
  
  message("Estimating intervals using LOD and Bayes method...")
  
  colnames_intervals <- c("phenotype", "pheno_num", "chr", "lod", 
                          "start.marker", "start.pos", 
                          "max.marker", "max.pos", 
                          "end.marker", "end.pos")
  
  # LOD method
  lodint <- matrix(ncol = 10, nrow = 0)
  colnames(lodint) <- colnames_intervals
  
  for (i in 1:dim(peaks)[1]){
    
    lodint.single <- lodint(out, chr = as.character(peaks$chr[i]), 
                            lodcolumn = as.numeric(peaks$pheno_num[i]), 
                            expandtomarkers = T, drop = 1)
    
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
  
  lodint <- distinct(lodint, pheno_num, chr, start.pos, max.pos, end.pos, .keep_all = T)
  
  # Find markers
  for (i in 1:dim(lodint)[1]){
    
    lodint[i,7] <- find.marker(cross, chr = lodint[i,3], pos = lodint[i,8])
    
  }
  
  # Bayes method
  bayesint <- matrix(ncol = 10, nrow = 0)
  colnames(bayesint) <- colnames_intervals
  
  for (i in 1:dim(peaks)[1]){
    
    bayesint.single <- bayesint(out, chr = as.character(peaks$chr[i]), 
                               lodcolumn = as.numeric(peaks$pheno_num[i]), 
                               expandtomarkers = T)
    
    bayesint.temp <- data.frame(phenotype = as.character(peaks$phenotype[i]),
                                pheno_num = as.numeric(peaks$pheno_num[i]), 
                                chr = as.character(peaks$chr[i]), 
                                lod = peaks$lod[i], 
                                start.marker = row.names(bayesint.single)[1], 
                                start.pos = bayesint.single[1,2], 
                                max.marker = row.names(bayesint.single)[2], 
                                max.pos = bayesint.single[2,2], 
                                end.marker = row.names(bayesint.single)[3], 
                                end.pos = bayesint.single[3,2])
    
    bayesint <- rbind(bayesint, bayesint.temp)
    
  }
  
  bayesint <- distinct(bayesint, pheno_num, chr, start.pos, max.pos, end.pos, .keep_all = T) 
  
  # Find markers
  for (i in 1:dim(lodint)[1]) {
    
    bayesint[i,7] <-  find.marker(cross, chr = bayesint[i,3], pos = bayesint[i,8])
    
  }
  
  
  
  # 7: QTL effects -------------------------------------------------------------
  
  message("Computing QTL effects")
  
  qtl.effects <- cbind(lodint[,c(1,3,4,7)], additive = NA, dominance = NA, ratio = NA)
  
  for (i in 1:dim(lodint)[1]){
    
    mar <- find.marker(cross, chr = lodint$chr[i], pos = lodint$max.pos[i])
    eff <- effectplot(cross, mname1 = mar, pheno.col = lodint$pheno_num[i] + 1, draw = F)
    qtl.effects$additive[i] <- (eff[["Means"]][1] - eff[["Means"]][4])/2
    qtl.effects$dominance[i] <-
      (eff[["Means"]][2] + eff[["Means"]][3])/2 - (eff[["Means"]][1] + eff[["Means"]][4])/2
    # effectplot(cross, mname1 = mar, pheno.col = lodint$pheno_num[i]+1)
    # plotPXG(cross, marker = mar, pheno.col = lodint$pheno_num[i]+1)
    
  }
  
  qtl.effects$ratio <- qtl.effects$dominance/abs(qtl.effects$additive) %>%
    mutate(heritability = 1 - 10^(-(2/as.numeric(nind))*lod))
  
  
  
  # 8: Saving data -------------------------------------------------------------

  # Creates the directories to save the files and the plots
  outputdir <- paste0(dir, prefixResults, "_Results/")
  outputplot <- paste0(dir, prefixResults, "_Plots/")
  dir.create(outputdir, showWarnings = F)
  dir.create(outputplot, showWarnings = F)
  
  message("Saving data...")
  
  # R Data
  setwd(dir)
  save.image(paste0(dir, prefixResults, ".RData"))
  
  # Save files
  write.table(lodint,
              file = paste(outputdir, prefixResults, ".qtl.LodIntervals.txt", sep = ""),
              sep = "\t", row.names = F, quote = F)
  write.table(bayesint,
              file = paste(outputdir, prefixResults, ".qtl.BayesIntervals.txt", sep = ""),
              sep = "\t", row.names = F, quote = F)
  write.table(out, file = paste(outputdir, prefixResults, ".lod.scores.txt", sep = ""),
              sep = "\t", row.names = T, quote = F)
  write.table(qtl.effects,
              file = paste(outputdir, prefixResults, ".qtl.effects.txt", sep = ""),
              sep = "\t", row.names = F, quote = F)
  
  
  
  # 9: Plots -------------------------------------------------------------------
  
  message(paste0("Saving plots at: ", outputdir, " and ", outputplot))
  
  source(paste0(dirfun, "02_QTL_SinglePlots.R"))
  genoPlot.pdf(cross)
  missGenoPlot.pdf(cross)
  plotPhenotypes.pdf(data.pheno)
  plotLOD.pdf(out)
  
  source(paste0(dirfun, "02_QTL_MapPlot.R") )
  plot.single.QTL(lodint, "Lod")
  plot.single.QTL(bayesint, "Bayes")
  
  message("Done!")
  
}



###### Example(s) ######
# Set arguments
# dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F1_Metabolomics/02_QTL_Analysis/CM8996/"
# dircross <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F1_Metabolomics/02_QTL_Analysis/CM8996/"
# dirfun <- "D:/OneDrive - CGIAR/GitHub/CassBio.QTL_analysis/R code/"
# locfile <- "CM8996_postmapping_v.2_2023.loc"
# mapfile <- "CM8996_LM_2023.map"
# phenofile <- "CM8996_metabolomic_QTL.qua"
# prefixResults <- "CM8996_QTL"
# ncores <- 50

# Run function
# single_qtl(dir, dircross, dirfun, locfile, mapfile, phenofile, prefixResults, ncores)
