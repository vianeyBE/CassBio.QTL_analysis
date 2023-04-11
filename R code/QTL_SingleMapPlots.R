# Short name: Single map plots
# Description: Creates PDF files with the chromosome, interval, and LOD profile
# Output: PDFs with the plots
# 
# Author: Vianey Barrera (vpbarrerae@gmail.com / v.barrera@cgiar.org)
# 
# Arguments:
#



plot.single.QTL <- function(qtldata, int_method){
  
  library(RColorBrewer)
  library(ggtext)
  
  options(warn = -1)
  
  # 1: Preparing information ---------------------------------------------------
  
  # Map data
  data.map <- matrix(ncol = 3, nrow = 0)
  colnames(data.map) <- c("mar", "pos", "chr")
  
  # Chromosome size
  chr.size <- matrix(ncol = 2, nrow = 0)
  colnames(chr.size) <- c("chromosome", "size")
  
  # LG names
  nameschro = names(cross[["geno"]])
  
  for (i in c(1:nchro)){
    
    # Map
    chr.map <- cross[["geno"]][[i]][["map"]][1,]
    temp.map <- data.frame(mar = names(chr.map), 
                          chr = (rep(nameschro[i], length(chr.map)[1])),
                          pos = chr.map, 
                          row.names = NULL)
    data.map <- rbind(data.map, temp.map)
    
    # Chromosome
    chr.size.temp <- data.frame(chromosome = nameschro[i], size = max(temp.map$pos))
    chr.size <- rbind(chr.size, chr.size.temp)
    
  }
  
  
  
  # 2: Reading each chromosome -------------------------------------------------

  for (LG in unique(qtldata$chr)){
    
    message(paste0("Creating plot for chromosome ", LG))
    
    data.map.qtl <- subset(data.map, chr == LG) # Get data of the map
    chr.size.qtl <- subset(chr.size, chromosome == LG ) # Get data of the size
    data.qtl <- subset(qtldata, chr == LG) # Get qtls on the LG
    mypal.qtl <- brewer.pal(n = dim(data.qtl)[1], name = "Dark2")[1:dim(data.qtl)[1]] 
    qtl.label <- as.data.frame(rbind(cbind(data.qtl$start.marker, data.qtl$start.pos),
                      cbind(data.qtl$max.marker, data.qtl$max.pos),
                      cbind(data.qtl$end.marker, data.qtl$end.pos)))
    colnames(qtl.label) <- c("markers", "pos")
    qtl.label$pos <- as.numeric(qtl.label$pos)
    
    

    # 3: Plot chromosome -------------------------------------------------------

    # Axis 
    wide.chr <- 0.75
    brk <- 5
    my_breaks <- seq(0, round(max(chr.size.qtl$size), digits = -1) + brk, by = brk)
    my_labs <- interleave(seq(0, round(max(chr.size.qtl$size), digits = -1), by = brk * 2), "")
    
    # A plot with the chromosome, SNPs, QTL intervals, and labels of flank SNPs inside the chromosome
    qtlmap.basic <- function(chr.size.qtl, data.map.qtl, data.qtl){
      
      qtl.map <- ggplot() + 
        scale_x_discrete(name = "", position = "top", expand = expansion(add = c(1, 1))) + # Axis
        scale_y_reverse(name = "Location (cM)", expand = expansion(mult = c(0, 0)),
                        breaks = my_breaks, labels = my_labs) +
        geom_rect(data = chr.size.qtl, 
                  aes(xmin = 1 - wide.chr / 3, xmax = 1 + wide.chr / 3,
                      ymin = 0, ymax = size + 2), 
                  alpha = 0, color = "grey20") + # Draw chromosomes
        geom_rect(data = data.map.qtl,
                  aes(xmin = c(rep(1, dim(data.map.qtl)[1])) - wide.chr / 3,
                      xmax = c(rep(1, dim(data.map.qtl)[1])) + wide.chr / 3,
                      ymin = pos - 0.001, ymax = pos + 0.001),
                  alpha = 0.65, color = "black", fill = "black") + # Draw SNPs
        geom_rect(data = data.qtl,
                  aes(xmin = c(rep(1,dim(data.qtl)[1])) - wide.chr / 3,
                      xmax = c(rep(1,dim(data.qtl)[1])) + wide.chr / 3,
                      ymin = start.pos, ymax = end.pos),
                  alpha = 0.2, fill = mypal.qtl, color = "black") + # Draw qtl intervals
        theme(panel.background = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.length = unit(.5, "cm"),
              axis.text.x = element_text(size = 24, face = "bold"),
              axis.text.y = element_text(size = 24, face = "bold", color = "black",
                                         margin = margin(0, 0, 0, 20)),
              text = element_text(size = 24, face = "bold"),
              axis.line.y.left = element_line(color = "black", size = 1)) +
        ggrepel::geom_text_repel(data = qtl.label, #Labels
                                 mapping = aes(x = c(rep(1, dim(qtl.label)[1])), 
                                             y = pos, label = markers),
                               size = 14, nudge_x = - 0.75, direction = "y")
      
      return(qtl.map)
    
    }
    
    qtl.map <- qtlmap.basic(chr.size.qtl, data.map.qtl, data.qtl)
    
    
    
    # 4: Plot interval and label -----------------------------------------------

    data.qtl$xvalue <- 1 + wide.chr * 0.4 + seq(from = 0, to = 0.1 * (dim(data.qtl)[1] - 1), by = 0.1)
    
    qtl_range <- list()
    qtl_text <- list()
    
    for (i in 1:dim(data.qtl)[1]){
      
      p <- geom_pointrange(data = data.qtl[i,],
                           aes(x = xvalue, y = max.pos, ymin = start.pos, ymax = end.pos),
                           size = 2, shape = 16, alpha = 1, fatten = 2,colour = mypal.qtl[i])
      q <- geom_richtext(data = data.qtl[i,],
                         aes(x = xvalue, y = end.pos, # it might be necessary to move label (!)
                             label = phenotype, angle = -90, hjust = -0.05),
                         label.color = NA, fill = NA, text.colour = mypal.qtl[i], size = 14)
      qtl_range[[i]] <- p
      qtl_text[[i]] <- q
    
    }
    
    qtl.map <- qtl.map + list(qtl_range) + list(qtl_text)
  
    
    
    # 5: LOD profile -----------------------------------------------------------
    
    qtl.lod.temp <- data.frame(loci =  row.names(out), chr = as.character(out$chr), pos = out[,2],
                               out[,c(data.qtl$phenotype)])
    
    qtl.lod.temp <- subset(qtl.lod.temp, chr == LG)
    
    plotLODmqm <- function(data.lod.temp){
      
      row.names(data.lod.temp <- c())
      maxlod <- max(as.data.frame(data.lod.temp[,4:dim(data.lod.temp)[2]]))
      data.lod <- melt(data.lod.temp, id.vars = c("loci", "chr", "pos"))
      data.lod.markers <- subset(data.lod, loci %in% data.map$mar)
        
      my_breaks <- seq(0, round(maxlod) + 0.5,by = 0.5)
      my_labs <- interleave(seq(0,round(maxlod), by = 1), "")
      # data.lod <- melt(data.lod.temp, id.vars = c("loci", "chr", "pos"))
        
        plot.lod <- ggplot(data.lod, aes(x = pos, y = value)) + 
          geom_line(size = 1 , aes(color = variable)) 
          geom_rug(data = data.lod.markers, aes(x = pos), sides = "b",
                   position = "jitter", alpha = 0.5, size = 0.25) + 
          scale_x_reverse(breaks = seq(0, round(max(chr.size.qtl$size), digits = -1) + 5, by = 10),
                          expansion(mult = c(0, 0.))) +
          scale_y_continuous(position = "right", breaks = my_breaks, lab = my_labs) +
          scale_color_manual(values = mypal.qtl) +
          labs(x = "", y = "") +
          coord_flip() +
          theme_classic() + 
          theme(strip.text = element_text(size = 15, face = "bold", family = font),
                strip.background = element_rect(colour = NA, fill = NA),
                strip.placement = "outside",
                axis.text = element_text(size = 25, face = "bold", family = font, 
                                     color = "black", margin = margin(0, 0, 0, 15)),
                text = element_text(size = 25,face = "bold", family = font),
                axis.line.y = element_line(color = "black", size = 1),
                legend.position = "bottom", legend.title = element_blank(),
                legend.box = "vertical", legend.text = element_text(size = 25),
                legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(20, 20, 20, 20))
    
        return(plot.lod)
          
    }
    
    plot.lod <- plotLODmqm(qtl.lod.temp)
    
    

    # 6: Combine and save plots ------------------------------------------------
    
    qtls.plot <- plot_grid(qtl.map, plot.lod, rel_widths = c(0.7, 1))
    
    save_plot(plot = qtls.plot, base_height = 30,
              filename = paste(outputplot, prefixResults, ".single_qtl.", LG, ".",
                               int_method, ".pdf", sep = ""))
    
  }
  
  message("QTL plots Done!")
  
}
