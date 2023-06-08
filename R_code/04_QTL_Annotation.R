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
# 1: Add a new pathway for common regions



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
    
    gff <- read.delim(gff, header = F, comment.char = "#") %>%
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
    
    gff <- read.delim(gff, header = F, comment.char = "#") %>%
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
    csv_L <- list()
    
    # Loop to find QTL files
    for (i in 1:length(names)){
      
      # Database handling 
      csv_L[[i]] <- read.delim(paste0(Wdir, "/", names[i])) %>%
        select(phenotype, chr, lod, start.marker, end.marker) %>%
        mutate(start = start.marker) %>%
        tidyr::separate(col = start, into = c("na", "Start"), sep = paste("_")) %>%
        mutate(end = end.marker) %>%
        tidyr::separate(col = end, into = c("na", "End"), sep = paste("_")) %>%
        select(phenotype, chr, lod, start.marker, end.marker, Start, End) %>%
        rename(Phenotype = phenotype, Chr = chr, LOD = lod,
               Start.Marker = start.marker, End.Marker = end.marker)
      
      # Progress bar
      cat('\r', i, ' files processed |', rep('=', i / 4), ifelse(i == length(names), '|\n', '>'), sep = '')
      
    }
    
    # Merge the data frames of the list in a single data frame
    s_QTL <- bind_rows(csv_L)
    s_QTL$Start <- as.numeric(s_QTL$Start)
    s_QTL$End <- as.numeric(s_QTL$End)
    
    # Generates the merged QTL data frame
    m_QTL <- s_QTL %>% group_by(Chr) %>%
      summarise("QTL.Start" = min(Start), "QTL.End" = max(End)) %>%
      mutate(QTL.Wide = QTL.End - QTL.Start)
    
  } else {
    
    # Read the files with the results
    setwd(Wdir)
    
    s_QTL <- read.csv(paste0(name)) %>%
      select(phenotype, chr, lod, start.marker, end.marker) %>%
      mutate(start = start.marker) %>%
      tidyr::separate(col = start, into = c("na", "Start"), sep = paste("_"), extra = "drop") %>%
      mutate(end = end.marker) %>%
      tidyr::separate(col = end, into = c("na", "End"), sep = paste("_"), extra = "drop") %>%
      select(phenotype, chr, lod, start.marker, end.marker, Start, End) %>%
      rename(Phenotype = phenotype, Chr = chr, LOD = lod,
             Start.Marker = start.marker, End.Marker = end.marker)
    
    s_QTL$Start <- as.numeric(s_QTL$Start)
    s_QTL$End <- as.numeric(s_QTL$End)
    
    m_QTL <- s_QTL %>% group_by(Chr) %>%
      summarise("QTL.Start" = min(Start), "QTL.End" = max(End)) %>%
      mutate(QTL.Wide = QTL.End - QTL.Start)
    
  }
  
  
  
  # 3: Match the QTL results with the gene annotation databases ----------------
  
  # Conditional: If there is at least one QTL to annotate, the function continues
  if (dim(s_QTL)[1] > 0){
    
    message(paste("There are", dim(s_QTL)[1], "QTLs to annotate"))
    
    
    
    # 3.1 Function pathway for single QTLs (s_QTL) -----------------------------
    
    # 3.1.1 Starts filtering the gff3 data -------------------------------------
    
    message("Obtaining gene information...")
    
    # Creates an empty list object to save the data
    s_gff_l <- list()
    
    # Filters and merges the gff3 data
    for (i in 1:nrow(s_QTL)){
      
      # Filter the gff3 file per chromosome
      data <- dplyr::filter(gff, s_QTL$Chr[i] == Chr)
      
      # Combine the different filtered data
      s_gff_l[[i]] <- rbind(
        
        # Looking for genes located inside the QTL
        dplyr::filter(data, s_QTL$Start[i] <= Start & s_QTL$Start[i] <= End &
                        s_QTL$End[i] >= Start & s_QTL$End[i] >= End) %>%
          mutate(Trait = s_QTL$Phenotype[i], LOD = s_QTL$LOD[i], Location = "Gene inside QTL",
                 QTL.Start = s_QTL$Start[i], QTL.End = s_QTL$End[i]) %>%
          rename(Gen.Start = Start, Gen.End = End),
        
        # Looking for genes that contain the QTL
        dplyr::filter(data, s_QTL$Start[i] >= Start & s_QTL$Start[i] <= End &
                        s_QTL$End[i] >= Start & s_QTL$End[i] <= End) %>%
          mutate(Trait = s_QTL$Phenotype[i], LOD = s_QTL$LOD[i], Location = "Gene inside QTL",
                 QTL.Start = s_QTL$Start[i], QTL.End = s_QTL$End[i]) %>%
          rename(Gen.Start = Start, Gen.End = End),
        
        # Looking genes down stream 
        dplyr::filter(data, s_QTL$Start[i] >= Start & s_QTL$Start[i] <= End) %>%
          mutate(Trait = s_QTL$Phenotype[i], LOD = s_QTL$LOD[i], Location = "Gene overlaps QTL",
                 QTL.Start = s_QTL$Start[i], QTL.End = s_QTL$End[i]) %>%
          rename(Gen.Start = Start, Gen.End = End),
        
        # Looking genes up stream
        dplyr::filter(data, s_QTL$End[i] >= Start & s_QTL$End[i] <= End) %>%
          mutate(Trait = s_QTL$Phenotype[i], LOD = s_QTL$LOD[i], Location = "Gene overlaps QTL",
                 QTL.Start = s_QTL$Start[i], QTL.End = s_QTL$End[i]) %>%
          rename(Gen.Start = Start, Gen.End = End)
        
      )
      
      # Progress bar
      cat('\r', i, ' files processed |', rep('=', i / 4), ifelse(i == nrow(s_QTL), '|\n', '>'), sep = '')
      
    }
    
    # Merge all the data frames of the list in a unique data frame
    s_gff_m <- bind_rows(s_gff_l) %>%
      select(Trait, Chr, LOD, QTL.Start, QTL.End, What, Name, Location, Gen.Start, Gen.End)
    
    
    
    # 3.2.2 Continues filtering the annotation data ----------------------------
    
    message("Obtaining gene annotation...")
    
    # Creates an empty list object to save the data
    s_annot_l <- list()
    
    # Filters and merges the annot data
    for (i in 1:nrow(s_gff_m)){
      
      # Filter the annotation data and creates new columns based on previous data
      s_annot_l[[i]] <- filter(annot,
                                Locus == s_gff_m$Name[i] |
                                  Trans == s_gff_m$Name[i] |
                                    Peptide == s_gff_m$Name[i]) %>%
        mutate(Trait =s_gff_m$Trait[i], Chr = s_gff_m$Chr[i], LOD = s_gff_m$LOD[i],
               QTL.Start = s_gff_m$QTL.Start[i], QTL.End = s_gff_m$QTL.End[i],
               What = s_gff_m$What[i], Location = s_gff_m$Location[i],
               Gen.Start = s_gff_m$Gen.Start[i], Gen.End = s_gff_m$Gen.End[i])
      
      # Progress bar
      cat('\r', i, ' rows processed |', rep('=', i / 10), ifelse(i == nrow(s_gff_m), '|\n', '>'), sep = '')
      
    }
    
    
    
    # 3.2.3 It ends by formatting the data frame -------------------------------
    
    # Merge the data frames of the list in a single data frame and modify it
    s_QTL_annotation <- bind_rows(s_annot_l) %>%
      select(Trait, Chr, LOD, QTL.Start, QTL.End, What, ID, Locus, Trans, Peptide, Location, 
             Gen.Start, Gen.End, GO, AT.name, AT.define)
    
    
    
    # 3.2 Function pathway for merged QTLs (m_QTL) -----------------------------
    
    # 3.2.1 Starts filtering the gff3 data -------------------------------------
    
    message("Obtaining gene information...")
    
    # Creates an empty list object to save the data
    m_gff_l <- list()
    
    # Filters and merges the gff3 data
    for (i in 1:nrow(m_QTL)){
      
      # Filter the gff3 file per chromosome
      data <- dplyr::filter(gff, m_QTL$Chr[i] == Chr)
      
      # Combine the different filtered data
      m_gff_l[[i]] <- rbind(
        
        # Looking for genes located inside the QTL
        dplyr::filter(data, m_QTL$QTL.Start[i] <= Start & m_QTL$QTL.Start[i] <= End &
                        m_QTL$QTL.End[i] >= Start & m_QTL$QTL.End[i] >= End) %>%
          mutate(QTL.Start = m_QTL$QTL.Start[i], QTL.End = m_QTL$QTL.End[i],
                 Location = "Gene inside QTL") %>%
          rename(Gen.Start = Start, Gen.End = End),
        
        # Looking for genes that contain the QTL
        dplyr::filter(data, m_QTL$QTL.Start[i] >= Start & m_QTL$QTL.Start[i] <= End &
                        m_QTL$QTL.End[i] >= Start & m_QTL$QTL.End[i] <= End) %>%
          mutate(QTL.Start = m_QTL$QTL.Start[i], QTL.End = m_QTL$QTL.End[i],
                 Location = "QTL inside gene") %>%
          rename(Gen.Start = Start, Gen.End = End),
        
        # Looking genes down stream 
        dplyr::filter(data, m_QTL$QTL.Start[i] >= Start & m_QTL$QTL.Start[i] <= End) %>%
          mutate(QTL.Start = m_QTL$QTL.Start[i], QTL.End = m_QTL$QTL.End[i],
                 Location = "Gene overlaps QTL") %>%
          rename(Gen.Start = Start, Gen.End = End),
        
        # Looking genes up stream
        dplyr::filter(data, m_QTL$QTL.End[i] >= Start & m_QTL$QTL.End[i] <= End) %>%
          mutate(QTL.Start = m_QTL$QTL.Start[i], QTL.End = m_QTL$QTL.End[i],
                 Location = "Gene overlaps QTL") %>%
          rename(Gen.Start = Start, Gen.End = End)
        
      )
      
      # Progress bar
      cat('\r', i, ' files processed |', rep('=', i / 1), ifelse(i == nrow(m_QTL), '|\n', '>'), sep = '')
      
    }
    
    # Merge the data frames of the list in a single data frame
    m_gff_m <- bind_rows(m_gff_l) %>%
      select(Chr, QTL.Start, QTL.End, What, Name, Location, Gen.Start, Gen.End)
    
    
    
    # 3.2.2 Continues filtering the annotation data ----------------------------
    
    message("Obtaining gene annotation...")
    
    # Creates an empty list object to save the data
    m_annot_l <- list()
    
    # Filters and merges the annot data
    for (i in 1:nrow(m_gff_m)){
      
      # Filter the annotation data and creates new columns based on previous data
      m_annot_l[[i]] <- filter(annot,
                                Locus == m_gff_m$Name[i] |
                                  Trans == m_gff_m$Name[i] |
                                    Peptide == m_gff_m$Name[i]) %>%
        mutate(Chr = m_gff_m$Chr[i], What = m_gff_m$What[i], Location = m_gff_m$Location[i],
               Gen.Start = m_gff_m$Gen.Start[i], Gen.End = m_gff_m$Gen.End[i],
               QTL.Start = m_gff_m$QTL.Start[i], QTL.End = m_gff_m$QTL.End[i])
      
      # Progress bar
      cat('\r', i, ' rows processed |', rep('=', i / 10), ifelse(i == nrow(m_annot_l), '|\n', '>'), sep = '')
      
    }
    
    
    
    # 3.2.3 It ends by formatting the data frame -------------------------------
    
    # Merge the data frames of the list in a single data frame and modify it
    m_QTL_annotation <- bind_rows(m_annot_l) %>% 
      select(Chr, QTL.Start, QTL.End, What, ID, Locus, Trans, Peptide, Location, 
             Gen.Start, Gen.End, GO, AT.name, AT.define)
    
    
    
    # 4: Save the outputs ------------------------------------------------------
    
    message("Saving output files as: 'QTL_annotation.csv' and 'QTL_merged_annotation.csv'")
    
    write.csv(s_QTL_annotation, file = "single_QTL_annotation.csv", quote = F, row.names = F)
    write.csv(m_QTL_annotation, file = "merged_QTL_annotation,csv", quote = F, row.names = F)
    
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
# gff <- "Mesculenta_305_v6.1/Mesculenta_305_v6.1.gene.gff3"
# version <- "6.1"
# recursive <- "F"

# Run function
# QTL_Annotation(Wdir, Ddir, name, wdyw, recursive)
