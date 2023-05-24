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
# pat: Enter the path of file names to look for. For example: QTL_LOD_Intervals.
# wdyw: Enter what are you looking for to annotate (Options: CDS, five_prime_UTR, gene, mRNA, three_prime_UTR).
# annot: Annotation details of the genes. txt file from the genome version used for alignment.
# GFF: gff3 file from the genome version used for alignment.
# version: You can choose between the genome of reference version 6.1 or 8.1 (Options: 6.1 or 8.1).



###### To do ######
# 1: Update single QTL annotation
# 2: Add QTL-Trait name



# 0: Function init -------------------------------------------------------------

QTL_Annotation <- function(Wdir, Ddir, pat, wdyw, annot, GFF, version){
  
  # 1: Load all the directories, info, and data --------------------------------
  
  # Load packages
  if (!require(tidyverse)) install.packages(tidyverse)
  library(tidyverse)
  
  # Load annotation data and re-organize it
  
  message("Loading required files to do the annotation")
  
  if (version == 6.1){
    
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
  
  
  
  # 2: Find all the CSVs with the results and filter them ----------------------
  
  message("Getting list of CSV files...")
  
  # Get the names of the files
  setwd(Wdir)
  names <- list.files(path = Wdir, pattern = paste0(pat, "."), all.files = F, full.names = F, recursive = T)
  
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
  QTL_s <- QTL %>% group_by(Chr) %>% summarise("QTL.Start" = min(Start), "QTL.End" = max(End))
  
  message(paste("There are", dim(QTL_s)[1], "QTLs after filtering"))
  
  
  # 3: Match the results with the gene annotation database ---------------------
  
  # Creates a conditional
  # If there is at least one result of the GWAS models, the annotation is done
  if (dim(QTL_s)[1] > 0){
    
    message("Retriving annotation...")
    
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
      annotL[[i]] <- rbind(
        
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
    annotLm <- bind_rows(annotL)
    
    message(paste("There are", dim(annotLm)[1], "genes that overlap with QTLs"))
    
    
    
    # 4: Formatting dataframe --------------------------------------------------
    
    message("Obtaining gene information...")
    
    # Creates the objects to perform the loop
    p <- 0
    gffL <- list()
    
    # Loop
    for (i in 1:nrow(annotLm)){
      
      # Iterator
      p <- p + 1
      
      #
      gffL[[i]] <- filter(annot,
                          Locus == annotLm$Name[p] |
                            Trans == annotLm$Name[p] |
                            Peptide == annotLm$Name[p]) %>%
        mutate(Chr = annotLm$Chr[p], 
               Gen.Start = annotLm$Gen.Start[p], Gen.End = annotLm$Gen.End[p],
               QTL.Start = annotLm$QTL.Start[p], QTL.End = annotLm$QTL.End[p], 
               Location = annotLm$Location[p])
      
      # Progress bar
      cat('\r', i, ' rows processed |', rep('=', i / 5), ifelse(i == nrow(annotLm), '|\n', '>'), sep = '')
      
    }
    
    # Merge the data frames of the list in a single data frame and modify it
    QTL_annotation <- bind_rows(gffL) %>% 
      select(Chr, Locus, Gen.Start, Gen.End, GO, AT.name, AT.define, QTL.Start, QTL.End, Location)
    
    
    
    # 5: Save output -----------------------------------------------------------
    
    message("Saving output file: 'QTL_Annotation.csv' in working directory")
    
    write.csv(QTL_annotation, file = "QTL_Annotation.csv", quote = F, row.names = F)
    
    message("Done! ")
    
  } else {
    
    message("No QTLs to annotate")
    
  }
  
}



###### Example(s) ######
# Set arguments
# Wdir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F1_Metabolomics/02_QTL_Analysis/CM8996/All/"
# Ddir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/00_Data/"
# pat <- "LodIntervals"
# wdyw <- "gene"
# annot <- "Mesculenta_305_v6.1/Mesculenta_305_v6.1.annotation_info.txt"
# GFF <- "Mesculenta_305_v6.1/Mesculenta_305_v6.1.gene.gff3"
# version <- "6.1"

# Run function
# QTL_Annotation(Wdir, Ddir, pat, wdyw)
