# Short name: Annotation for QTL results
# Description (1): Get the annotation of the gen containing the SNP using the rqtl/QTL maping results
# Description (2): Retrieves the closest genes down and up stream
# Output: Basename of the output file
#
# Authors: Camilo E. Sanchez (c.e.sanchez@cgiar.org) and Vianey Barrera-Enriquez (vpbarrera@gmail.com)
#
# Arguments:
# Wdir: Name of the directory that contains the GAPIT results. For example: home/user/folder.
# Ddir: Directory where is located the annotation files (annot, GFF files).
# name: Enter the path or the name of file names to look for. For example: QTL_LOD_Intervals.
# wdyw: Enter what are you looking for to annotate (Options: CDS, five_prime_UTR, gene, mRNA, three_prime_UTR).
# annot: Annotation details of the genes. txt file from the genome version used for alignment.
# GFF: gff3 file from the genome version used for alignment.
# version: You can choose between the genome of reference version 6.1 or 8.1 (Options: 6.1 or 8.1).
# recursive: A Boolean string that determines if the function should perform a recursive search or not (Default = F).



###### To do ######
# 1: Include trait name to QTL annotation (QTL object)



# 0: Function init -------------------------------------------------------------

QTL_Annotation <- function(Wdir, Ddir, name, wdyw, annot, GFF, version, recursive = F){
  
  
  
  # 1: Load all the directories, info, and data --------------------------------
  
  # Load packages
  if (!require(tidyverse)) install.packages(tidyverse)
  library(tidyverse)
  
  # Load annotation data and re-organize it
  
  if (version == 6.1){
    
    message("Loading required files to do the annotation")
    
    # Set working directory where files are located
    setwd(Ddir)
    
    annot <- read.delim(annot, header = F) %>%
      rename(ID = 1, Locus = 2, Trans = 3, Peptide = 4, GO = 10, AT.name = 12, AT.define = 13) %>%
      select(ID, Locus, Trans, Peptide, GO, AT.name, AT.define)
    
    GFF <- read.delim(GFF, header = F, comment.char = "#") %>%
      rename(Chr = V1, What = V3, Start = V4, End = V5) %>%
      tidyr::separate(col = V9, into = c("ID", "na"), sep = ";") %>%
      tidyr::separate(col = ID, into = c("na2", "na3"), sep = "=") %>%
      tidyr::separate(col = na3, into = c("Name", "na4"), sep = ".v") %>%
      select(Chr, What, Start, End, Name) %>%
      mutate(Chr = recode(Chr, Chromosome01 = 1, Chromosome02 = 2, Chromosome03 = 3, Chromosome04 = 4,
                          Chromosome05 = 5, Chromosome06 =  6, Chromosome07 = 7, Chromosome08 = 8, 
                          Chromosome09 = 9, Chromosome10 = 10, Chromosome11 = 11, Chromosome12 = 12,
                          Chromosome13 = 13, Chromosome14 = 14, Chromosome15 = 15, Chromosome16 = 16, 
                          Chromosome17 = 17, Chromosome18 = 18)) %>%
      na.omit() %>% 
      dplyr::filter(What %in% wdyw)
    
  } else {
    
    message("Loading required files to do the annotation")
    
    # Set working directory where files are located
    setwd(Ddir)
    
    annot <- read.delim(annot, header = T) %>%
      rename(ID = 1, Locus = 2, Trans = 3, Peptide = 4, GO = 10, AT.name = 11, AT.define = 12) %>%
      select(ID, Locus, Trans, Peptide, GO, AT.name, AT.define)
    
    GFF <- read.delim(GFF, header = F, comment.char = "#") %>%
      rename(Chr = V1, What = V3, Start = V4, End = V5) %>%
      tidyr::separate(col = V9, into = c("ID", "na"), sep = ";") %>%
      tidyr::separate(col = ID, into = c("na2", "na3"), sep = "=") %>%
      tidyr::separate(col = na3, into = c("Name", "na4"), sep = ".v") %>%
      select(Chr, What, Start, End, Name) %>%
      mutate(Chr = recode(Chr, Chromosome01 = 1, Chromosome02 = 2, Chromosome03 = 3, Chromosome04 = 4,
                          Chromosome05 = 5, Chromosome06 =  6, Chromosome07 = 7, Chromosome08 = 8, 
                          Chromosome09 = 9, Chromosome10 = 10, Chromosome11 = 11, Chromosome12 = 12,
                          Chromosome13 = 13, Chromosome14 = 14, Chromosome15 = 15, Chromosome16 = 16, 
                          Chromosome17 = 17, Chromosome18 = 18)) %>%
      na.omit() %>% 
      dplyr::filter(What %in% wdyw)
    
  }
  
  
  
  # 2: Find the CSV(s) with the results and filter them ------------------------
  
  if (recursive == T){
    
    message("Getting list of CSV files...")
    
    # Get the names of the files
    setwd(Wdir)
    names <- list.files(path = Wdir, pattern = paste0(name, "."), all.files = F, full.names = F, recursive = T)
    
    message(paste("QTL mapping files found:", length(names)))
    message("Reading QTL mapping files...")
    
    # Create an empty list and the variable to iterate
    # Read, rename, modify and save all the csv files in the list
    csvL <- list()
    p <- 0
    
    # Loop to find QTL files
    for (i in 1:length(names)){
      
      # Iterator
      p <- p + 1
      
      # Database handling 
      csvL[[i]] <- read.delim(paste0(Wdir, "/", names[p])) %>%
        select(phenotype, chr, lod, start.marker, end.marker) %>%
        mutate(start = start.marker) %>%
        tidyr::separate(col = start, into = c("na", "Start"), sep = paste("_")) %>%
        mutate(end = end.marker) %>%
        tidyr::separate(col = end, into = c("na", "End"), sep = paste("_")) %>%
        select(phenotype, chr, lod, start.marker, end.marker, Start, End) %>%
        rename(Phenotype = phenotype, Chr = chr, LOD = lod, Start.Marker = start.marker, End.Marker = end.marker)
      
      # Progress bar
      cat('\r', i, ' files processed |', rep('=', i / 4), ifelse(i == length(names), '|\n', '>'), sep = '')
      
    }
    
    # Merge the data frames of the list in a single data frame
    QTL <- bind_rows(csvL)
    QTL$Start <- as.numeric(QTL$Start)
    QTL$End <- as.numeric(QTL$End)
    
    QTL_s <- QTL %>% group_by(Chr) %>%
      summarise("QTL.Start" = min(Start), "QTL.End" = max(End)) %>%
      mutate(QTL.Wide = QTL.End - QTL.Start)
    
  } else {
    
    # Read the files with the results
    setwd(Wdir)
    
    QTL <- read.csv(paste0(name)) %>%
      select(phenotype, chr, lod, start.marker, end.marker) %>%
      mutate(start = start.marker) %>%
      tidyr::separate(col = start, into = c("na", "Start"), sep = paste("_"), extra = "drop") %>%
      mutate(end = end.marker) %>%
      tidyr::separate(col = end, into = c("na", "End"), sep = paste("_"), extra = "drop") %>%
      select(phenotype, chr, lod, start.marker, end.marker, Start, End) %>%
      rename(Phenotype = phenotype, Chr = chr, LOD = lod, Start.Marker = start.marker, End.Marker = end.marker)
    
    QTL$Start <- as.numeric(QTL$Start)
    QTL$End <- as.numeric(QTL$End)
    
    QTL_s <- QTL %>% group_by(Chr) %>%
      summarise("QTL.Start" = min(Start), "QTL.End" = max(End)) %>%
      mutate(QTL.Wide = QTL.End - QTL.Start)
    
  }
  
  
  
  # 3: Match the QTL results with the gene annotation database -----------------
  
  # Creates a conditional
  # If there is at least one result of the GWAS models, the annotation continues
  if (dim(QTL)[1] > 0){
    
    message(paste("There are", dim(QTL)[1], "QTLs to annotate"))
    
    
    
    # 3.1 Function path for individual QTLs (QTL) ------------------------------
    
    # 3.1.1  ----------------------------
    
    # Creates the objects to perform the loop
    p <- 0
    annotL <- list()
    
    # Loop
    for (i in 1:nrow(QTL)){
      
      # Iterator
      p <- p + 1
      
      # 
      data <- dplyr::filter(GFF, QTL$Chr[p] == Chr)
      
      # 
      annotL[[i]] <- rbind(
        
        # Looking for genes located inside the QTL
        dplyr::filter(data, QTL$Start[p] <= Start & QTL$Start[p] <= End &
                        QTL$End[p] >= Start & QTL$End[p] >= End) %>%
          mutate(QTL.Start = QTL$Start[p], QTL.End = QTL$End[p],
                 Location = "Gene inside QTL") %>%
          rename(Gen.Start = Start, Gen.End = End),
        
        # Looking for genes that contain the QTL
        dplyr::filter(data, QTL$Start[p] >= Start & QTL$Start[p] <= End &
                        QTL$End[p] >= Start & QTL$End[p] <= End) %>%
          mutate(QTL.Start = QTL$Start[p], QTL.End = QTL$End[p],
                 Location = "QTL inside gene") %>%
          rename(Gen.Start = Start, Gen.End = End),
        
        # Looking genes down stream 
        dplyr::filter(data, QTL$Start[p] >= Start & QTL$Start[p] <= End) %>%
          mutate(QTL.Start = QTL$Start[p], QTL.End = QTL$End[p],
                 Location = "Gene overlaps QTL") %>%
          rename(Gen.Start = Start, Gen.End = End),
        
        # Looking genes up stream
        dplyr::filter(data, QTL$End[p] >= Start & QTL$End[p] <= End) %>%
          mutate(QTL.Start = QTL$Start[p], QTL.End = QTL$End[p],
                 Location = "Gene overlaps QTL") %>%
          rename(Gen.Start = Start, Gen.End = End)
        
      )
      
      # Progress bar
      cat('\r', i, ' files processed |', rep('=', i / 4), 
          ifelse(i == nrow(QTL), '|\n', '>'), sep = '')
      
    }
    
    # Merge the data frames of the list in a single data frame and modify it
    annotLs <- bind_rows(annotL)
    
    message(paste("There are", dim(annotLs)[1], "genes that overlap with QTLs"))
    
    
    
    # 3.2.2: Formatting dataframe ----------------------------------------------
    
    message("Obtaining gene information...")
    
    # Creates the objects to perform the loop
    p <- 0
    gffL <- list()
    
    # Loop
    for (i in 1:nrow(annotLs)){
      
      # Iterator
      p <- p + 1
      
      #
      gffL[[i]] <- filter(annot,
                          Locus == annotLs$Name[p] |
                            Trans == annotLs$Name[p] |
                              Peptide == annotLs$Name[p]) %>%
        mutate(Chr = annotLs$Chr[p], 
               Gen.Start = annotLs$Gen.Start[p], Gen.End = annotLs$Gen.End[p],
               QTL.Start = annotLs$QTL.Start[p], QTL.End = annotLs$QTL.End[p], 
               Location = annotLs$Location[p])
      
      # Progress bar
      cat('\r', i, ' rows processed |', rep('=', i / 10), ifelse(i == nrow(annotLs), '|\n', '>'), sep = '')
      
    }
    
    # Merge the data frames of the list in a single data frame and modify it
    QTL_a <- bind_rows(gffL) %>% 
      select(Chr, Locus, Gen.Start, Gen.End, GO, AT.name, AT.define, QTL.Start, QTL.End, Location)
    
    
    
    # 3.2 Function path for merged QTLs (QTL_s) --------------------------------
    
    # 3.2.1  ------------------------------
    
    # Creates the objects to perform the loop
    p <- 0
    annotL <- list()
    
    # Loop
    for (i in 1:nrow(QTL_s)){
      
      # Iterator
      p <- p + 1
      
      # 
      data <- dplyr::filter(GFF, QTL_s$Chr[p] == Chr)
      
      # 
      annotLm[[i]] <- rbind(
        
        # Looking for genes located inside the QTL
        dplyr::filter(data, QTL_s$QTL.Start[p] <= Start & QTL_s$QTL.Start[p] <= End &
                        QTL_s$QTL.End[p] >= Start & QTL_s$QTL.End[p] >= End) %>%
          mutate(QTL.Start = QTL_s$QTL.Start[p], QTL.End = QTL_s$QTL.End[p],
                 Location = "Gene inside QTL") %>%
          rename(Gen.Start = Start, Gen.End = End),
        
        # Looking for genes that contain the QTL
        dplyr::filter(data, QTL_s$QTL.Start[p] >= Start & QTL_s$QTL.Start[p] <= End &
                        QTL_s$QTL.End[p] >= Start & QTL_s$QTL.End[p] <= End) %>%
          mutate(QTL.Start = QTL_s$QTL.Start[p], QTL.End = QTL_s$QTL.End[p],
                 Location = "QTL inside gene") %>%
          rename(Gen.Start = Start, Gen.End = End),
        
        # Looking genes down stream 
        dplyr::filter(data, QTL_s$QTL.Start[p] >= Start & QTL_s$QTL.Start[p] <= End) %>%
          mutate(QTL.Start = QTL_s$QTL.Start[p], QTL.End = QTL_s$QTL.End[p],
                 Location = "Gene overlaps QTL") %>%
          rename(Gen.Start = Start, Gen.End = End),
        
        # Looking genes up stream
        dplyr::filter(data, QTL_s$QTL.End[p] >= Start & QTL_s$QTL.End[p] <= End) %>%
          mutate(QTL.Start = QTL_s$QTL.Start[p], QTL.End = QTL_s$QTL.End[p],
                 Location = "Gene overlaps QTL") %>%
          rename(Gen.Start = Start, Gen.End = End)
        
      )
      
      # Progress bar
      cat('\r', i, ' files processed |', rep('=', i / 4), ifelse(i == nrow(QTL_s), '|\n', '>'), sep = '')
      
    }
    
    # Merge the data frames of the list in a single data frame
    annotLmm <- bind_rows(annotLm)
    
    message(paste("There are", dim(annotLm)[1], "genes that overlap with QTLs"))
    
    
    
    # 3.2.2: Formatting dataframe ----------------------------------------------
    
    message("Obtaining gene information...")
    
    # Creates the objects to perform the loop
    p <- 0
    gffL <- list()
    
    # Loop
    for (i in 1:nrow(annotLmm)){
      
      # Iterator
      p <- p + 1
      
      #
      gffL[[i]] <- filter(annot,
                          Locus == annotLmm$Name[p] |
                            Trans == annotLmm$Name[p] |
                            Peptide == annotLmm$Name[p]) %>%
        mutate(Chr = annotLmm$Chr[p], 
               Gen.Start = annotLmm$Gen.Start[p], Gen.End = annotLmm$Gen.End[p],
               QTL.Start = annotLmm$QTL.Start[p], QTL.End = annotLmm$QTL.End[p], 
               Location = annotLmm$Location[p])
      
      # Progress bar
      cat('\r', i, ' rows processed |', rep('=', i / 5), ifelse(i == nrow(annotLm), '|\n', '>'), sep = '')
      
    }
    
    # Merge the data frames of the list in a single data frame and modify it
    QTL_m_a <- bind_rows(gffL) %>% 
      select(Chr, Locus, Gen.Start, Gen.End, GO, AT.name, AT.define, QTL.Start, QTL.End, Location)
    
    
    
    # 4: Save the outputs --------------------------------------------------------
    
    message("Saving output files as: 'QTL_annotation.csv' and 'QTL_merged_annotation.csv'")
    
    write.csv(QTL_a, file = "QTL_merged_annotation.csv", quote = F, row.names = F)
    write.csv(QTL_m_a, file = "QTL_merged_annotation.csv", quote = F, row.names = F)
    
    message("Done!")
    
  } else {
    
    message("No QTLs to annotate")
    
  }
  
}



###### Example(s) ######
# Set arguments
# Wdir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F1_Metabolomics/02_QTL_Analysis/CM8996/All_metabolites/"
# Wdir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F1_Metabolomics/02_QTL_Analysis/CM8996/Significant_ones/"
# Ddir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/00_Data/"
# name <- "LodIntervals"
# name <- "QTL_results_heritability.csv"
# wdyw <- "gene"
# annot <- "Mesculenta_305_v6.1/Mesculenta_305_v6.1.annotation_info.txt"
# GFF <- "Mesculenta_305_v6.1/Mesculenta_305_v6.1.gene.gff3"
# version <- "6.1"
# recursive <- "F"

# Run function
# QTL_Annotation(Wdir, Ddir, name, wdyw, recursive)
