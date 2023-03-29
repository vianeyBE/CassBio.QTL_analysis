library("viridis")
library("reshape2")
library("ggplot2")
library("lemon")
library("colorspace")
library("cowplot") 


#### Plots for Rqtl #####

font = "sans"
interleave <- function(x,y){
  lx <- length(x)
  ly <- length(y)
  n <- max(lx,ly)
  as.vector(rbind(rep(x, length.out=n), rep(y, length.out=n)))
}

# PALETTE #
mypal.gen = plasma(ngeno)

# GENOTYPE #
genoPlot.pdf = function(cross){
  genotype.plot <- genoPlot(cross)
  file = paste(outputplot, prefixResults, ".genotype-plot.pdf", sep = "")
  pdf(file=file, paper='usr', width=10, height=8)
  print(genotype.plot)
  dev.off()
}
genoPlot = function(cross){
  # Data 
  data = matrix(ncol = 4, nrow=0)
  colnames(data) = c("Individual", "marker", "genotype", "chromosome")
  for (i in c(1:nchro)){
    data.chr = cross[["geno"]][[i]][["data"]]
    temp = melt(data.chr)
    temp = cbind(temp, chromosome = c(rep(i, dim(temp)[1])))
    colnames(temp) = c("Individual", "marker", "genotype", "chromosome")
    data = rbind(data, temp)
  }
  
  # Plot
  genotype.plot =
    ggplot(data, aes(x=marker, y=as.factor(Individual))) + 
      geom_raster(aes(fill=as.factor(genotype)), alpha = 0.8)+
      facet_grid(~chromosome, scales = "free", space = "free") +
      theme(legend.position="none",
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.text  = element_blank(),
            text = element_text(face = "bold", family=font),
            strip.text = element_text( face = "bold", family=font),
            strip.background = element_rect(colour = NA, fill = NA),
            panel.spacing = unit(0, "lines"),
            panel.background = element_rect(fill = "white", colour = "grey20")
      ) + 
      xlab("Makers") + ylab("Individuals") + labs(fill = "Genotypes") +
      scale_fill_manual(values = mypal.gen,
                        name ="", 
                        #labels = c("a", "b","c"),
                        na.value="grey30"
                        )
  return(genotype.plot)
}

# MISSING GENOTYPES #
missGenoPlot.pdf = function(cross){
  missing.plot <- missGenoPlot(cross)
  file = paste(outputplot, prefixResults, ".missing-genotype-plot.pdf", sep = "")
  pdf(file=file, paper='usr', width=10, height=8)
  print(missing.plot)
  dev.off()
}
missGenoPlot = function(cross){
  data = matrix(ncol = 4, nrow=0)
  colnames(data) = c("Individual", "marker", "genotype", "chromosome")
  for (i in c(1:nchro)){
    data.chr = cross[["geno"]][[i]][["data"]]
    temp = melt(data.chr)
    temp = cbind(temp, chromosome = c(rep(i, dim(temp)[1])))
    colnames(temp) = c("Individual", "marker", "genotype", "chromosome")
    data = rbind(data, temp)
  }
  
  missing.plot = 
    ggplot(data, aes(x=marker, y=as.factor(Individual))) + 
    geom_raster(aes(fill=as.factor(genotype)))+
    facet_grid(~chromosome, scales = "free", space = "free") +
    theme(legend.position="none",
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text  = element_blank(),
          text = element_text(face = "bold",family=font),
          strip.text = element_text(face = "bold", family=font),
          strip.background = element_rect(colour = NA, fill = NA),
          panel.spacing = unit(0, "lines"),
          panel.background = element_rect(fill = "white", colour = "black")
    ) + 
    xlab("Makers") + ylab("Individuals") + labs(fill = "Genotypes") +
    scale_fill_manual(values = c(rep("White", ngeno)),
                      #name ="", labels = c("a", "b","c"),
                      na.value="black")
  return(missing.plot)
}

# PHENOTYPES #
plotPhenotypes.pdf <- function(data.pheno) {
  
  plot_list <- list()
  
  pdf(paste(outputplot, prefixResults, ".phenotype-all-2.pdf", sep = ""), onefile=TRUE)
  for (i in 2:length(phenotypes)){
    
    #tmp <- 
      print(
      ggplot(data.pheno, aes(x=data.pheno[,i])) + 
      # histogram and density
        geom_histogram(aes( y=..density..),
                       colour = "black", fill = "white") +
        geom_density(alpha=.2, fill="#FF6666", size = 1, color = "#FF6666") + 
        geom_rug(aes(x = data.pheno[,i], y = 0), position = position_jitter(height = 0)) +
      # Aesthetic
        theme(panel.background = element_blank(),
              axis.text.x = element_text(face = "bold", family = font, color = "black"),
              axis.text.y = element_text(face = "bold", family = font, color = "black", margin=margin(0,0,0,20)),
              text = element_text(face = "bold", family = font),
              axis.line.y.left = element_line(color = "black", size = 1),
              axis.line.x.bottom = element_line(color = "black", size = 1),
              panel.border=element_blank(),
              axis.line = element_line()
        ) + 
      # Labs and limits
        xlab(phenotypes[i]) + ylab("") +
        xlim(min(na.omit(data.pheno[,i]))*0.95, max(na.omit(data.pheno[,i]))*1.05 ) +
        ylim(0, max(density(na.omit(data.pheno[,i]))$y)*1.20) +
        coord_capped_cart(left='both')
      )
    
    #plot_list[[i-1]] <- print(tmp)
  }
  
  #phenotypes.all <- plot_grid(plotlist=plot_list, ncol=3)
  
  #pdf(paste(outputplot, prefixResults, ".phenotype-all-2.pdf", sep = ""), onefile=TRUE)
  #print(plot_grid(plotlist=plot_list))
  dev.off()
  
  #phenotypes.all <- plot_grid(plotlist=plot_list)

  #save_plot(plot=phenotypes.all, filename = paste(outputplot, prefixResults, ".phenotype-all.pdf", sep = ""),base_height = 30)
}
plotPhenotypes =  function(data.pheno){
  for (i in 2:length(phenotypes)){
    print(
      ggplot(data.pheno, aes_string(x = data.pheno[,i])) + 
        geom_histogram(aes(y=..density..),  
                       colour = "black", fill = "white") +
        geom_density(alpha=.2, fill="#FF6666", size = 1, color = "#FF6666") + 
        geom_rug(aes_string(x = data.pheno[,i], y = 0), position = position_jitter(height = 0)) +
        theme(panel.background = element_blank(),
              axis.text.x = element_text(face = "bold", family = font, color = "black"),
              axis.text.y = element_text(face = "bold", family = font, color = "black", margin=margin(0,0,0,20)),
              text = element_text(face = "bold", family = font),
              axis.line.y.left = element_line(color = "black", size = 1),
              axis.line.x.bottom = element_line(color = "black", size = 1),
              panel.border=element_blank(),
              axis.line = element_line()
        ) + xlab(phenotypes[i]) + ylab("") +
        coord_capped_cart(left='both')
    )
  }
}

# LOD SCORE PROFILE #
plotLODCore = function(data.lod.temp){
  #out$chr = factor(out$chr, levels = c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", 
  #                                             "LG9", "LG10", "LG11", "LG12",  "LG13",  "LG14",  "LG15", 
  #                                             "LG16", "LG17","LG18"))
  #data.lod.temp = cbind(loci = row.names(out), data.frame(out))
  row.names(data.lod.temp) = c()
  maxlod = max(as.data.frame(data.lod.temp[,4:dim(data.lod.temp)[2]]))
  data.lod = melt(data.lod.temp, id.vars = c("loci", "chr", "pos"))
  data.lod.markers = subset(data.lod, loci %in% data.map$mar)
  
  my_breaks = seq(0, round(maxlod)+0.5,by=0.5)
  my_labs = interleave(seq(0,round(maxlod), by=1), "")
  data.lod = melt(data.lod.temp, id.vars = c("loci", "chr", "pos"))
  
  plot.lod= 
    ggplot(data.lod, aes(x=pos, y=value)) + 
    geom_line(size = 0.5 , aes(color = variable))+
    facet_grid(~chr,  scales = "free",
               space = "free", switch = "x" 
    ) +  
    geom_rug(data=data.lod.markers, 
             aes(x=pos), 
             sides = "b",
             position = "jitter",
             alpha = 0.5,
             size = 0.25
    ) + 
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0)),
                       breaks = my_breaks,
                       lab = my_labs,
                       sec.axis = dup_axis(breaks=my_breaks)) +
    coord_capped_cart(bottom='both') +
    geom_hline(yintercept = maxlod*1.05, size = 1, color = "black") +
    geom_hline(yintercept = 0, size = 1) +
    theme_classic() + 
    theme(#panel.background = element_blank(),
      panel.spacing = unit(0, "lines"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.ticks.y.left = element_line(size=1),
      axis.text.y.right = element_blank(),
      strip.text = element_text(face = "bold", family=font),
      strip.background = element_rect(colour = NA, fill = NA),
      strip.placement = "outside",
      #axis.ticks.length.y.left = unit(.1, "cm"),
      #axis.text.x = element_text(size=20, face = "bold", family = font, color = "black"),
      axis.text.y = element_text(face = "bold", family = font, color = "black", margin=margin(0,0,0,15)),
      text = element_text(face = "bold", family = font),
      axis.line.y = element_line(color = "black", size = 1),
      legend.position ="bottom",
      legend.title = element_blank(),
      legend.box="vertical", legend.margin=margin()
      #axis.line.x.bottom = element_line(color = "black", size = 1),
      #panel.border=element_blank(),
      #axis.line = element_line()
    ) + xlab("") + ylab("LOD score") + 
    scale_color_manual(values =  divergingx_hcl(n = npheno-1, palette = "Roma", alpha = 1)) +
    guides(col = guide_legend(ncol=3,byrow=TRUE))
  return(plot.lod)
}

plotLOD =  function(out){ ## Plot all phenotypes
  
  levelsLG <-  c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", 
               "LG9", "LG10", "LG11", "LG12",  "LG13",  "LG14",  "LG15", 
               "LG16", "LG17","LG18")
  levelschr <- c(1:18)
  out$chr <- factor(out$chr, levels = levelschr)
  
  data.lod.temp = cbind(loci = row.names(out), data.frame(out)) 
  plot.lod <- plotLODCore(data.lod.temp)
  return(plot.lod)
}
plotLOD.pdf = function(out){
  plot.lod <- plotLOD(out)
  file = paste(outputplot, prefixResults, ".lod-profile-all.pdf", sep = "")
  pdf(file=file, width=25, height=15)
  print(plot.lod)
  dev.off()
}

plotLOD.one = function(out, phenotypeName){ #Plot an specific  phenotype
  out$chr = factor(out$chr, levels = c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", 
                                       "LG9", "LG10", "LG11", "LG12",  "LG13",  "LG14",  "LG15", 
                                       "LG16", "LG17","LG18"))
  data.lod.temp = cbind(loci = row.names(out), data.frame(out))
  data.lod.temp = data.lod.temp[, c("loci", "chr", "pos", phenotypeName)]
  plot.lod <- plotLODCore(data.lod.temp)
  return(plot.lod)
}
plotLOD.one.pdf = function(out, phenotypeName){
  plot.lod <- plotLOD.one(out, phenotypeName)
  file = paste(outputplot, prefixResults, ".lod-profile.", phenotypeName, ".pdf", sep = "")
  pdf(file=file, width=25, height=15)
  print(plot.lod)
  dev.off()
  
}
