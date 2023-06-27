# Short name: QTL boxplots
# Description: Plot a boxplot for a SNP: x is genotype, y is phenotype. It is possible to add extra information about the samples using the labelfile
# Output: A single PDF file containing the boxplot of the SNPs
#
# Authors: Camilo E. Sánchez-Sarria (c.e.sanchez@cgiar.org) and Vianey Barrera-Enriquez (vpbarrerae@gmail.com)
#
# Arguments:
# prefix: Name that will have the output
# dir: Directory where is located the data and where output will be save
# phenofile: Phenotype data in csv format. Columns: samples, first one named as: Taxa'.
# genofile: Genotype data in hapmap format
# snp_list: CSV file with three columns:
#       01: SNPS. List of SNPS to plot, name should be the same as in the geno data
#       02: trait. Name of the trait as in the pheno data
#       03: xlabel. Name of the trait to be included as label
# recursive: 
# labelfile: CSV file with two columns (optional)
#        01: taxa. Sample names
#        02: label. Label or category to add at the plot
# order: 



##### To do #####
# 1: 



# 0: Function init -------------------------------------------------------------

QTL_Boxplot <- function(prefix, dir, phenofile, genofile, code, snp_list, recursive, labelfile, order){

  
  
  # 1: Load all the directories, info, and data --------------------------------
  
  # Load packages
  if (!require(tidyverse)) install.packages(tidyverse)
  if (!require(janitor)) install.packages(janitor)
  if (!require(Biostrings)) install.packages(Biostrings)
  if (!require(hrbrthemes)) install.packages(hrbrthemes)
  if (!require(ggsignif)) install.packages(ggsignif)
  if (!require(RColorBrewer)) install.packages(RColorBrewer)
  if (!require(openxlsx)) install.packages(openxlsx)
  
  library(tidyverse)
  library(janitor)
  library(Biostrings)
  library(hrbrthemes)
  library(ggsignif)
  library(RColorBrewer)
  library(openxlsx)
  
  # Set working directory
  setwd(dir)
  
  # Load files
  pheno <- read.csv(phenofile, head = T)
  geno <- read.delim(genofile, head = F)
  geno <- geno %>% janitor::row_to_names(1) %>%
    tidyr::separate(col = 1, into = c("chr", "posit"), sep = "_")
  
  # Conditional for how is coded the SNPs in the hapmap
  if (code == "S"){
    
    # Data frame modification
    geno <- geno %>% mutate(chrs = str_replace(chr, "S", "C")) %>%
      unite(SNPS, chrs, posit, sep = "_") %>%
      select(-c(1)) %>%
      remove_rownames() %>%
      column_to_rownames(var = 'SNPS')
    
  } else {
    
    if (code == "Chr"){
      
      # Data frame modification
      geno <- geno  %>% mutate(chrs = str_replace(chr, "Chromosome", "C")) %>%
        unite(SNPS, chrs, posit, sep = "_") %>%
        select(-c(1)) %>%
        remove_rownames() %>%
        column_to_rownames(var = 'SNPS')
      
    } 
  }

  # Conditional to read the SNPs list
  if (recursive == T){
    
    # Get the names of the files
    message("Getting list of txt files...")
    names <- list.files(path = dir, pattern = snp_list, all.files = F, full.names = F, recursive = T)
    
    # Informative messages
    message(paste("QTL LOD intervals files found:", length(names)))
    message("QTL LOD intervals files...")
    
    # Create an empty list and the variable to iterate
    # Read, rename, modify and save all the csv files in the list
    csv_L <- list()
    
    # Loop to find QTL files
    for (i in 1:length(names)){
      
      # Database handling 
      csv_L[[i]] <- read.delim(paste0(dir, names[i])) %>%
        select(max.marker, phenotype) %>%
        tidyr::separate(col = max.marker, into = c("chr", "pos", "nu"), sep = "_") %>%
        unite(SNPS, chr, pos, sep = "_") %>%
        mutate(xlabel = phenotype, trait = phenotype) %>%
        select(SNPS, trait, xlabel)
      
      # Progress bar
      cat('\r', i, ' files processed |', rep('=', i / 2), ifelse(i == length(names), '|\n', '>'), sep = '')
      
    }
    
    # Creates the single QTL data frame (s_QTL)
    # Merge the data frames of the list in a single data frame
    snp_list <- bind_rows(csv_L)
    
    # Sort the SNPs alphabetically
    snp_list$SNPS <- sort(snp_list$SNPS)
    
    } else {
      
      if (recursive == F) {
      
        snp_list <- "CM8996_metabolomic_results_plots.csv"
        
      # Read the SNPs list from a csv file
      snp_list <- read.csv(snp_list, header = T) %>%
        select(start.marker, max.marker, end.marker, phenotype)
      
      # Creates three data frames based on start, max, and end marker
      snp1 <- snp_list %>% select(start.marker, phenotype) %>% 
        dplyr::rename(marker = start.marker)
      snp2 <- snp_list %>% select(max.marker, phenotype) %>% 
        dplyr::rename(marker = max.marker)
      snp3 <- snp_list %>% select(end.marker, phenotype) %>% 
        dplyr::rename(marker = end.marker)
      
      # Merge data bases
      snp_list <- rbind(snp1, snp2, snp3)
      
      # Deletes the nucleotide info
      snp_list <- snp_list %>% 
        tidyr::separate(col = marker, into = c("chr", "pos", "nu"),sep = "_") %>%
        unite(SNPS, chr, pos, sep = "_") %>%
        mutate(xlabel = phenotype, trait = phenotype) %>%
        select(SNPS, trait, xlabel) %>%
        distinct(SNPS, .keep_all = T)
      
      # Sort the SNPs alphabetically
      snp_list$SNPS <- sort(snp_list$SNPS)
      rm(snp1, snp2, snp3)
      
    }
    
  }
  
  # Informative messages
  message("SNPs to Plot: ", length(unique(snp_list$SNPS)))
  
  # Creates a new nucleotide codification
  iupac <- Biostrings::IUPAC_CODE_MAP
  iupac["A"]<- "AA"
  iupac["C"]<- "CC"
  iupac["G"]<- "GG"
  iupac["T"]<- "TT"
  iupac["V"]<- NA
  iupac["H"]<- NA
  iupac["D"]<- NA
  iupac["B"]<- NA
  iupac["N"]<- NA
  
  
  
  # 2: Starting the process of drawing -----------------------------------------
  
  if (!is.null(labelfile)){
    
    # 2.1: With labels ---------------------------------------------------------
    
    message("Labels were provided, they will be included in the boxplot")
    
    
    
    # 2.1.1: Prepare the data --------------------------------------------------
    # Load label file and modify it
    label <- read.csv(labelfile, header = T,
                      col.names = c('Taxa', 'Levels'),
                      colClasses = c('character', 'factor'))
    
    message("Labels to colour data: ", paste(levels(label$Levels), collapse = ', '))
    
    # Creates a color palette
    palette <- brewer.pal(n = length(levels(label$Levels)), name = "Dark2")
    
    # Creates a copy of label file
    matrix_GT <- label
    
    # Create a new workbook
    message("Creating workbook")
    wb <- createWorkbook()
    
    # Creates the PDF file where boxplots will paste
    pdf(paste0(prefix, ".SNPs_boxplot.pdf", sep = ''), onefile = T)
    
    # Creates the boxplot for each trait
    for (i in 1:dim(snp_list)[1]){
      
      # Obtain each one of the snp, trait, and label
      snpname <- snp_list$SNPS[i]
      trait <- snp_list$trait[i]
      x_label <- snp_list$xlabel[i]
      
      # Join label and genotype data
      matrix_GT <- merge(label, t(geno[snpname, -c(1:10)]), by.x = 'Taxa', by.y = 'row.names')
      
      # Join label-genotype and pheno data
      data <- merge(matrix_GT, pheno[,c("Taxa", paste0("X", trait))], by = 'Taxa')
      
      # Replace IUPAC to nucleotides
      data$snp <- as.factor(iupac[data[, snpname]])
      
      # Delete NAs
      data <- na.omit(data)
      
      # Make a list of comparisons
      if (length(unique(data$snp)) == 2){
        
        comp <- as.list(data.frame(matrix(combn(unique(data$snp), 2), ncol = 1)))
        
      } else{
        
        comp <- as.list(data.frame(matrix(combn(unique(data$snp), 2), ncol = 3)))
        
      }
      
      # Re order the levels
      if (!is.null(order)){
        
        data$Levels <- factor(data$Levels, levels = order)
        
      }
      
      # Add the data to the workbook
      addWorksheet(wb, sheetName = paste0(colnames(data)[3]))
      writeData(wb, sheet = paste0(colnames(data)[3]), data, colNames = T)
      
      # 2.1.2: Draws the boxplots ----------------------------------------------
      
      plot <- ggplot(data, aes(x = snp, y = get(paste0("X", trait)))) +
        geom_violin(fill = "gray80", color = "white", width = 0.5, alpha = 0.5) +
        geom_boxplot(fill = "gray80", alpha = 0.75, width = 0.1) +
        geom_jitter(aes(color = Levels), size = 3, height = 0, width = 0.2, alpha = 1) +
        scale_color_manual(values = c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")) +
        geom_signif(comparisons = comp, map_signif_level = T) +
        labs(x = paste0("SNP: ", snpname), y = paste0("Trait: ", x_label)) +
        theme_classic() +
        theme(legend.position = "bottom", legend.title = element_blank(),
              legend.text =  element_text(size = 12, color = "black"),
              axis.text = element_text(size = 12, color = "black"),
              axis.title = element_text(size = 12, color = "black"),
              axis.title.x = element_text(margin = margin(t = 10)),
              axis.title.y = element_text(margin = margin(r = 10)))
        
      print(plot)
      
      # Progress bar
      cat('\r', i, ' SNPs processed |', rep('=', i/2), ifelse(i == dim(snp_list)[1], '|\n',  '>'), sep = '')
      
    }
    
    dev.off()
    
    # Save outputs
    saveWorkbook(wb, paste0(prefix, ".trait_x_snp.xlsx"), overwrite = T)
    
    dev.off()
    
    # Informative messages
    message(paste0("Outputs:\n\n", 
                   "1: ", prefix, ".trait_x_snp.xlsx\n",
                   "2: ", prefix, ".SNPs_boxplot.pdf\n\n",
                   "Were stored in: \n\n", dir))
    
  } else {
    
    # 2.2: Without labels ------------------------------------------------------
    
    message("No labels provided")
    
    # Creates the PDF file where boxplots will paste
    pdf(paste0(prefix, ".SNPs_boxplot.pdf", sep = ''), onefile = T)
    
    # Creates the boxplot for each trait
    for (i in 1:dim(snps)[1]){
      
      # Obtain each one of the snp, trait, and label
      snpname <- snp_list$SNPS[i]
      trait <- snp_list$trait[i]
      x_label <- snp_list$xlabel[i]
      
      # 2.2.1: Prepare the data ------------------------------------------------
      # Join pheno and geno data
      data <- merge(pheno, t(geno[snpname, -c(2:11)]), by.x = 'Taxa', by.y = 'row.names')
      
      # Replace IUPAC to nucleotides
      data$snp <- as.factor(iupac[data[, snpname]])
      
      # Delete the NAs
      data <- na.omit(data)
      
      # Make a list of comparisons
      if (length(unique(data$snp)) == 2){
        
        comp <- as.list(data.frame(matrix(combn(unique(data$snp), 2), ncol = 1)))
        
      } else{
        
        comp <- as.list(data.frame(matrix(combn(unique(data$snp), 2), ncol = 3)))
        
      }
      
      # Number of levels to plot
      n <- length(levels(data$snp))
      
      
      
      # 2.2.2: Draws the boxplots ----------------------------------------------
      
      plot <- ggplot(data, aes(x = snp, y = get(paste0("X", trait)))) +
        geom_violin(fill = "gray80", color = "gray80", width = 0.5, alpha = 0.5) +
        geom_boxplot(fill = c("#1b9e77", "#d95f02", "#7570b3")[1:length(levels(data$snp))],
                     alpha = 0.75, width = 0.2) +
        geom_jitter(size = 2, height = 0, width = 0.1, alpha = 1) +
        geom_signif(comparisons = comp, map_signif_level = T) +
        labs(x = paste0("SNP: ", snpname), y = paste0("Trait: ", x_label)) +
        theme_classic() +
        theme(legend.position = "none", 
              axis.text = element_text(size = 12, color = "black"),
              axis.title = element_text(size = 12, color = "black"),
              axis.title.x = element_text(margin = margin(t = 10)),
              axis.title.y = element_text(margin = margin(r = 10))) 
      
      print(plot)
      
      # Progress bar
      cat('\r', i, ' SNPs processed |', rep('=', i/2), ifelse(i == dim(snps)[1], '|\n',  '>'), sep = '')
      
    }
    
    dev.off()
    
  }
  
  
  
  # 3: Function ends -----------------------------------------------------------
  message("Done!")
  
}



###### Example(s) ######
# Set arguments
 prefix <- "F1_CM8996"
 dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F1_Metabolomics/02_QTL_Analysis/CM8996/"
 phenofile <- "CM8996_metabolomic.csv"
 genofile <- "CM8996.final.hmp.txt"
 code <- "Chr"
 labelfile <- "CM8996_labels.csv"
 order <- c("S", "IS", "I", "IR", "R")
 
 # If recursive == False
 snp_list <- "CM8996_metabolomic_results_plots.csv"
 recursive <- F
 
 # # If recursive == True
 snp_list <- ".LodIntervals"
 recursive <- T



# Run function
# QTL_Boxplot(prefix, dir, phenofile, genofile, code, snp_list, recursive, labelfile, order)
