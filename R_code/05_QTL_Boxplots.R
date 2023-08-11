# Short name: QTL boxplots
# Description: Plot a boxplot for a SNP: x is genotype, y is phenotype.
#              It is possible to add extra information about the samples using the labelfile.
# Output: A single PDF file containing the boxplot of the SNPs
#
# Authors: Camilo E. Sánchez-Sarria (c.e.sanchez@cgiar.org) and Vianey Barrera-Enriquez (vpbarrerae@gmail.com).
#
# Arguments:
# prefix: Name that will have the output.
# dir: Directory where is located the data and where output will be save.
# phenofile: Phenotype data in csv format. Columns: samples. Rows: individuals.
# genofile: Genotype data in hapmap format.
# snpList: CSV file with three columns:
#      01: List of SNPS to plot, name should be the same as in the geno data.
#      02: Name of the trait as in the pheno data.
# recursive: Boolean that indicates if there is to perform a recursive search of files (Options: 'T' or 'F').
# labelfile: CSV file with two columns (optional).
#        01: taxa. Sample names.
#        02: label. Label or category to add at the plot.
# order: A vector string indicating the order to plot given the labelfile (optional).



##### To do #####
# Everything is good!



# 0: Function init -------------------------------------------------------------

QTL_Boxplot <- function(prefix, dir, phenofile, genofile, snpList, recursive, labelfile = NULL, 
                        order = NULL){
  
  
  
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
  
  # Informative message
  message("Reading pheno and geno files")
  
  # Load files
  pheno <- read.csv(phenofile, head = T)
  colnames(pheno)[1] <- "Taxa"
  
  geno <- read.delim(genofile, head = F)
  geno <- geno %>% janitor::row_to_names(1) %>%
    tidyr::separate(col = 1, into = c("chr", "posit"), sep = "_")
  
  # Re-code the SNPs in the hapmap to 'C'
  if (grepl("S", geno[["chr"]])[1] == T){
    
    # Informative message
    message("SNPs are coded as 'S'\n\n", "Recoding them to 'C'")
    
    # Data frame modification
    geno <- geno %>% mutate(chrs = str_replace(chr, "S", "C")) %>%
      unite(SNPS, chrs, posit, sep = "_") %>%
      select(-c(1)) %>%
      remove_rownames() %>%
      column_to_rownames(var = 'SNPS')
    
  } else {
    
    if (grepl("Chr", geno[["chr"]])[1] == T){
      
      # Informative message
      message("SNPs are coded as 'Chr'\n\n", "Recoding them to 'C'")
      
      # Data frame modification
      geno <- geno  %>% mutate(chrs = str_replace(chr, "Chromosome", "C")) %>%
        unite(SNPS, chrs, posit, sep = "_") %>%
        select(-c(1)) %>%
        remove_rownames() %>%
        column_to_rownames(var = 'SNPS')
      
    } else {
      
      message("SNPs are coded correctly as 'C'")
      
    }
    
  }
  
  # Conditional to read the SNPs list
  if (recursive == T){
    
    # Get the names of the files
    message("Getting list of QTL LOD intervals txt files...\n\n")
    names <- list.files(path = dir, pattern = snpList, all.files = F, full.names = F, recursive = T)
    
    # Informative messages
    message(paste0("QTL LOD intervals files found: ", length(names)),
            "\n\nReading QTL LOD intervals files...")
    
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
        mutate(trait = phenotype) %>%
        select(SNPS, trait)
      
      # Progress bar
      cat('\r', i, ' files processed |', rep('=', i / 2), ifelse(i == length(names), '|\n', '>'), sep = '')
      
    }
    
    # Merge the data frames of the list in a single data frame
    snpList <- bind_rows(csv_L)
    
    # Order the snpList by chromosome
    snpList <- snpList[order(snpList$SNPS), ]
    
    } else {
      
      if (recursive == F){
      
      # Read the SNPs list from a csv file
      snpList <- read.csv(snpList, header = T) %>%
        select(start.marker, max.marker, end.marker, phenotype)
      
      # Creates three data frames based on start, max, and end marker
      snp1 <- snpList %>% select(start.marker, phenotype) %>% 
        dplyr::rename(marker = start.marker)
      snp2 <- snpList %>% select(max.marker, phenotype) %>% 
        dplyr::rename(marker = max.marker)
      snp3 <- snpList %>% select(end.marker, phenotype) %>% 
        dplyr::rename(marker = end.marker)
      
      # Merge data bases
      snpList <- rbind(snp1, snp2, snp3)
      
      # Deletes the nucleotide info
      snpList <- snpList %>% 
        tidyr::separate(col = marker, into = c("chr", "pos", "nu"), sep = "_") %>%
        unite(SNPS, chr, pos, sep = "_") %>%
        mutate(trait = phenotype) %>%
        select(SNPS, trait)
      
      # Order the snpList by chromosome
      snpList <- snpList[order(snpList$SNPS), ]
      
    }
    
  }
  
  # Informative messages
  message("SNPs to plot: ", length(unique(snpList$SNPS)), "\n\n",
          "Traits to plot:\n", paste(unique(snpList$trait), collapse = '\n'))
  
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
    
    # Creates a copy of label file
    matrix_GT <- label
    
    # Creates an empty list to storage the data
    data_snp <- list()
    
    # Creates the PDF file where boxplots will paste
    pdf(paste0(prefix, ".SNPs_boxplot.pdf", sep = ''), onefile = T)
    
    # Creates the boxplot for each trait
    for (i in 1:dim(snpList)[1]){
      
      # Obtain each one of the snp, trait, and label
      snpname <- snpList$SNPS[i]
      trait <- snpList$trait[i]
      
      # Join label and genotype data
      matrix_GT <- merge(label, t(geno[snpname, -c(1:10)]), by.x = 'Taxa', by.y = 'row.names')
      
      # Join label-genotype and pheno data
      data <- merge(matrix_GT, pheno[, c("Taxa", paste0("X", trait))], by = 'Taxa')
      
      # Replace IUPAC to nucleotides
      data$snp <- as.factor(iupac[data[, snpname]])
      
      # Delete NAs
      data <- na.omit(data)
      
      # Change column names
      colnames(data)[3] <- "na"
      colnames(data)[4] <- "Values"
      colnames(data)[5] <- "SNP"
      
      # Database handling
      data <- data %>% mutate(Trait = trait, Name = snpname) %>%
        select(Taxa, Levels, Values, SNP, Trait, Name)
      
      # Save the data in the list
      data_snp[[i]] <- data
      
      # Make a list of comparisons
      if (length(unique(data$SNP)) == 2){
        
        comp <- as.list(data.frame(matrix(combn(unique(data$SNP), 2), ncol = 1)))
        
      } else {
        
        if (length(unique(data$SNP)) == 1){
          
          message(paste0("For the trait: ", trait, " and the SNP: ", snpname,
                         " there are not possible comparasions to plot"))
          
        } else {
          
          comp <- as.list(data.frame(matrix(combn(unique(data$SNP), 2), ncol = 3)))
          
        }
        
      }
      
      # Re order the levels
      if (!is.null(order)){
        
        data$Levels <- factor(data$Levels, levels = order)
        
      }
      
      
      
      # 2.1.2: Draws the boxplots ----------------------------------------------
      
      # Plot depending of the amount of haplotypes
      if (length(unique(data$SNP)) == 1){
        
        # Plot for unique haplotypes
        plot <- ggplot(dataj, aes(SNP, Values)) +
          geom_violin(fill = "gray80", color = "white", width = 0.5, alpha = 0.5) +
          geom_boxplot(fill = "gray80", alpha = 0.75, width = 0.1) +
          geom_jitter(aes(color = Levels), size = 2, height = 0, width = 0.2, alpha = 1) +
          scale_color_brewer(palette = "Spectral") +
          labs(x = paste0("SNP: ", snpname), y = paste0("Trait: ", trait)) +
          theme_classic() +
          theme(legend.position = "bottom", legend.title = element_blank(),
                legend.text =  element_text(size = 12, color = "black"),
                axis.text = element_text(size = 12, color = "black"),
                axis.title = element_text(size = 12, color = "black"),
                axis.title.x = element_text(margin = margin(t = 10)),
                axis.title.y = element_text(margin = margin(r = 10)))
        
      } else {
        
        # Plot for several haplotypes
        plot <- ggplot(dataj, aes(SNP, Values)) +
          geom_violin(fill = "gray80", color = "white", width = 0.5, alpha = 0.5) +
          geom_boxplot(fill = "gray80", alpha = 0.75, width = 0.1) +
          geom_jitter(aes(color = Levels), size = 2, height = 0, width = 0.2, alpha = 1) +
          scale_color_brewer(palette = "Spectral") +
          geom_signif(comparisons = comp, map_signif_level = T) +
          labs(x = paste0("SNP: ", snpname), y = paste0("Trait: ", trait)) +
          theme_classic() +
          theme(legend.position = "bottom", legend.title = element_blank(),
                legend.text =  element_text(size = 12, color = "black"),
                axis.text = element_text(size = 12, color = "black"),
                axis.title = element_text(size = 12, color = "black"),
                axis.title.x = element_text(margin = margin(t = 10)),
                axis.title.y = element_text(margin = margin(r = 10)))
        
      }
      
      print(plot)
      
      # Progress bar
      cat('\r', i, ' SNPs processed |', rep('=', i/2), ifelse(i == dim(snpList)[1], '|\n',  '>'), sep = '')
      
    }
    
    dev.off()
    
    # Merge the list in one dataframe
    data_snp_M <- bind_rows(data_snp)
    
    # Create an excel workbook to save the data
    message("Creating excel workbook to store the data")
    wb <- createWorkbook()
    addWorksheet(wb, sheetName = paste0(prefix))
    writeData(wb, sheet = paste0(prefix), data_snp_M, colNames = T)
    saveWorkbook(wb, paste0(prefix, ".trait_x_snp.xlsx"), overwrite = T)
    
    # Informative messages
    message(paste0("Outputs:\n\n", 
                   "1: ", prefix, ".trait_x_snp.xlsx\n",
                   "2: ", prefix, ".SNPs_boxplot.pdf\n\n",
                   "Were stored in: \n\n", dir))
    
  } else {
    
    # 2.2: Without labels ------------------------------------------------------
    
    # Informative message
    message("No labels provided")
    
    # Creates an empty list to storage the data
    data_snp <- list()
    
    # Creates the PDF file where boxplots will paste
    pdf(paste0(prefix, ".SNPs_boxplot.pdf", sep = ''), onefile = T)
    
    # Creates the boxplot for each trait
    for (i in 1:dim(snpList)[1]){
      
      
      
      # 2.2.1: Prepare the data ------------------------------------------------
      
      # Obtain each one of the snp, trait, and label
      snpname <- snpList$SNPS[i]
      trait <- snpList$trait[i]
      
      # Join pheno and geno data
      data <- merge(pheno[, c("Taxa", paste0("X", trait))], t(geno[snpname, -c(2:11)]), 
                    by.x = "Taxa", by.y = "row.names")
      
      # Replace IUPAC to nucleotides
      data$snp <- as.factor(iupac[data[, snpname]])
      
      # Delete the NAs
      data <- na.omit(data)
      
      # Change column names
      colnames(data)[2] <- "Values"
      colnames(data)[3] <- "na"
      colnames(data)[4] <- "SNP"
      
      # Database handling
      data <- data %>% mutate(Trait = trait, Name = snpname) %>%
        select(Taxa, Values, SNP, Trait, Name)
      
      # Make a list of comparisons
      if (length(unique(data$SNP)) == 2){
        
        comp <- as.list(data.frame(matrix(combn(unique(data$SNP), 2), ncol = 1)))
        
      } else {
        
        if (length(unique(data$SNP)) == 1){
          
          message(paste0("For the trait: ", trait, " and the SNP: ", snpname,
                         " there are not possible comparasions to plot"))
          
        } else {
          
          comp <- as.list(data.frame(matrix(combn(unique(data$SNP), 2), ncol = 3)))
          
        }
        
      }
      
      # Save the data in the list
      data_snp[[i]] <- data
      
      
      
      # 2.2.2: Draws the boxplots ----------------------------------------------
      
      # Plot depending of the amount of haplotypes
      if (length(unique(data$SNP)) == 1){
        
        # Plot for unique haplotypes
        plot <- ggplot(data, aes(SNP, Values)) +
          geom_violin(fill = "gray80", color = "white", width = 0.5, alpha = 0.5) +
          geom_boxplot(aes(fill = SNP), alpha = 0.75, width = 0.1) +
          geom_jitter(size = 2, height = 0, width = 0.2, alpha = 1) +
          geom_signif(comparisons = comp, map_signif_level = T) +
          labs(x = paste0("SNP: ", snpname), y = paste0("Trait: ", trait)) +
          theme_classic() +
          theme(legend.position = "none", 
                axis.text = element_text(size = 12, color = "black"),
                axis.title = element_text(size = 12, color = "black"),
                axis.title.x = element_text(margin = margin(t = 10)),
                axis.title.y = element_text(margin = margin(r = 10)))
        
      } else {
        
        # Plot for several haplotypes
        plot <- ggplot(data, aes(SNP, Values)) +
          geom_violin(fill = "gray80", color = "white", width = 0.5, alpha = 0.5) +
          geom_boxplot(aes(fill = SNP), alpha = 0.75, width = 0.1) +
          geom_jitter(size = 2, height = 0, width = 0.2, alpha = 1) +
          scale_fill_brewer(palette = "Dark2") +
          geom_signif(comparisons = comp, map_signif_level = T) +
          labs(x = paste0("SNP: ", snpname), y = paste0("Trait: ", trait)) +
          theme_classic() +
          theme(legend.position = "none", 
                axis.text = element_text(size = 12, color = "black"),
                axis.title = element_text(size = 12, color = "black"),
                axis.title.x = element_text(margin = margin(t = 10)),
                axis.title.y = element_text(margin = margin(r = 10)))
        
      }
      
      print(plot)
      
      # Progress bar
      cat('\r', i, ' SNPs processed |', rep('=', i/2), ifelse(i == dim(snpList)[1], '|\n',  '>'), sep = '')
      
    }
    
    dev.off()
    
    # Merge the list in one dataframe
    data_snp_M <- bind_rows(data_snp)
    
    # Create an excel workbook to save the data
    message("Creating excel workbook to store the data")
    wb <- createWorkbook()
    addWorksheet(wb, sheetName = paste0(prefix))
    writeData(wb, sheet = paste0(prefix), data_snp_M, colNames = T)
    saveWorkbook(wb, paste0(prefix, ".trait_x_snp.xlsx"), overwrite = T)
    
    # Informative messages
    message(paste0("Outputs:\n\n",
                   "1: ", prefix, ".trait_x_snp.xlsx\n",
                   "2: ", prefix, ".SNPs_boxplot.pdf\n\n",
                   "Were stored in: \n\n", dir))
    
  }
  
  
  
  # 3: Function ends -----------------------------------------------------------
  
  message("Done!")
  
}



# Example(s) -------------------------------------------------------------------
# Set arguments
 prefix <- "F1_CM8996"
 dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F1_Metabolomics/02_QTL_Analysis/CM8996/"
 phenofile <- "CM8996_metabolomic.csv"
 genofile <- "CM8996.final.hmp.txt"
# labelfile <- "CM8996_labels.csv"
# order <- c("S", "IS", "I", "IR", "R")
 
# If recursive == False
# snpList <- "CM8996_metabolomic_results_plots.csv"
# recursive <- F
 
# If recursive == True
# snpList <- ".LodIntervals"
# recursive <- T



# Run function -----------------------------------------------------------------
# QTL_Boxplot(prefix, dir, phenofile, genofile, snpList, recursive, labelfile = NULL, order = NULL)
