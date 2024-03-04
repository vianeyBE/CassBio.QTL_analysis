# Short name: Annotation for QTL results
# Description (1): Get the annotation of the gen containing the SNP using the rqtl/QTL maping results
# Description (2): Retrieves the closest genes down and up stream
# Output: Basename of the output file
#
# Authors: Camilo E. Sanchez (c.e.sanchez@cgiar.org) and Vianey Barrera-Enriquez (vpbarrera@gmail.com)
#
# Arguments:
# Ddir: Directory where is located the annotation files (annot, GFF files).
# annot: Annotation details of the genes. txt file from the genome version used for alignment.
# gff: gff3 file from the genome version used for alignment.
# version: You can choose between the genome of reference version 6.1 or 8.1 (Options: 6.1 or 8.1).
# wdyw: Enter what are you looking for to annotate (Options: CDS, five_prime_UTR, gene, mRNA, three_prime_UTR).
# prefix: Prefix name that indicates how the outputs will be saved.
# Wdir: Name of the directory that contains the GAPIT results. For example: home/user/folder.
# name: Enter the path or the name of file names to look for. For example: QTL_LOD_Intervals.
# recursive: A Boolean string that determines if the function should perform a recursive search or not (Default = F).



###### To do ######
# 1: Modify location text
# 2: Add trait name in c_QTL
# 3: Modify version selection



# 0: Function init -------------------------------------------------------------

QTL_Annotation <- function(Ddir, annot, gff, version, wdyw, prefix, Wdir, name, recursive){
  
  
  
  # 1: Load all the directories, info, and data --------------------------------
  
  # Load packages
  if (!require(tidyverse)) install.packages(tidyverse)
  if (!require(openxlsx)) install.packages(openxlsx)
  
  library(tidyverse)
  library(openxlsx)
  
  # Load annotation data and re-organize it
  if (version == 6.1){
    
    # Informative message
    message("Loading required files to do the annotation")
    
    # Read and modify the files
    annot <- read.delim(paste0(Ddir, annot), header = F) %>%
      rename(ID = 1, Locus = 2, Trans = 3, Peptide = 4, GO = 10, AT.name = 12, AT.define = 13) %>%
      select(ID, Locus, Trans, Peptide, GO, AT.name, AT.define)
    
    gff <- read.delim(paste0(Ddir, gff), header = F, comment.char = "#") %>%
      rename(Chr = V1, What = V3, Start = V4, End = V5) %>%
      tidyr::separate(col = V9, into = c("ID", "na"), sep = ";", extra = "drop") %>%
      tidyr::separate(col = ID, into = c("na2", "na3"), sep = "=", extra = "drop") %>%
      tidyr::separate(col = na3, into = c("Name", "na4"), sep = ".v", extra = "drop") %>%
      select(Chr, What, Start, End, Name) %>%
      dplyr::filter(grepl("Chromosome", Chr)) %>%
      mutate(Chr = recode(Chr, Chromosome01 = 1, Chromosome02 = 2, Chromosome03 = 3, Chromosome04 = 4,
                          Chromosome05 = 5, Chromosome06 =  6, Chromosome07 = 7, Chromosome08 = 8, 
                          Chromosome09 = 9, Chromosome10 = 10, Chromosome11 = 11, Chromosome12 = 12,
                          Chromosome13 = 13, Chromosome14 = 14, Chromosome15 = 15, Chromosome16 = 16, 
                          Chromosome17 = 17, Chromosome18 = 18)) %>%
      dplyr::filter(What %in% wdyw)
    
    # Informative message
    message("Annotation files loaded succesfully")
    
  } else {
    
    # Informative message
    message("Loading required files to do the annotation")

    # Read and modify the files
    annot <- read.delim(paste0(Ddir, annot), header = T) %>%
      rename(ID = 1, Locus = 2, Trans = 3, Peptide = 4, GO = 10, AT.name = 11, AT.define = 12) %>%
      select(ID, Locus, Trans, Peptide, GO, AT.name, AT.define)
    
    gff <- read.delim(paste0(Ddir, gff), header = F, comment.char = "#") %>%
      rename(Chr = V1, What = V3, Start = V4, End = V5) %>%
      tidyr::separate(col = V9, into = c("ID", "na"), sep = ";", extra = "drop") %>%
      tidyr::separate(col = ID, into = c("na2", "na3"), sep = "=", extra = "drop") %>%
      tidyr::separate(col = na3, into = c("Name", "na4"), sep = ".v", extra = "drop") %>%
      select(Chr, What, Start, End, Name) %>%
      dplyr::filter(grepl("Chromosome", Chr)) %>%
      mutate(Chr = recode(Chr, Chromosome01 = 1, Chromosome02 = 2, Chromosome03 = 3, Chromosome04 = 4,
                          Chromosome05 = 5, Chromosome06 =  6, Chromosome07 = 7, Chromosome08 = 8, 
                          Chromosome09 = 9, Chromosome10 = 10, Chromosome11 = 11, Chromosome12 = 12,
                          Chromosome13 = 13, Chromosome14 = 14, Chromosome15 = 15, Chromosome16 = 16, 
                          Chromosome17 = 17, Chromosome18 = 18)) %>%
      dplyr::filter(What %in% wdyw)
    
    # Informative message
    message("Annotation files loaded succesfully")
    
  }
  
  
  
  # 2: Find the CSV(s) with the results and creates the data frames ------------
  if (recursive == T){
    
    # Set working directory
    setwd(Wdir)
    
    # Get the names of the files
    message("Getting list of CSV files...")
    names <- list.files(path = Wdir, pattern = paste0(name, "."), all.files = F, full.names = F, recursive = T)
    
    # Informative messages
    message(paste("QTL mapping files found:", length(names)))
    message("Reading QTL mapping files...")
    
    # Create an empty list to store data
    csv_L <- list()
    
    # Loop to find, read, rename, modify and save all QTL files
    for (i in 1:length(names)){
      
      # Database handling 
      csv_L[[i]] <- read.delim(paste0(Wdir, names[i])) %>%
        select(phenotype, chr, lod, start.marker, end.marker) %>%
        mutate(start = start.marker) %>%
        tidyr::separate(col = start, into = c("na", "Start"), sep = paste("_"), extra = "drop") %>%
        mutate(end = end.marker) %>%
        tidyr::separate(col = end, into = c("na", "End"), sep = paste("_"), extra = "drop") %>%
        select(phenotype, chr, lod, start.marker, end.marker, Start, End) %>%
        rename(Phenotype = phenotype, Chr = chr, LOD = lod,
               Start.Marker = start.marker, End.Marker = end.marker)
      
      # Progress bar
      cat('\r', i, ' files processed |', rep('=', i / 4), ifelse(i == length(names), '|\n', '>'), sep = '')
      
    }
    
    # Creates the single QTL data frame (s_QTL)
    # Merge the data frames of the list in a single data frame
    s_QTL <- bind_rows(csv_L)
    s_QTL$Start <- as.numeric(s_QTL$Start)
    s_QTL$End <- as.numeric(s_QTL$End)
    
    
    
    # Creates the converged QTL data frame (c_QTL)
    # Creates am empty list to save the QTLs converge regions per chromosomes 
    best <- list()
    
    # Find the regions in which converge most of QTLs per chromosome
    for (p in unique(s_QTL$Chr)){
      
      # Filter the QTL data per chromosome
      temp <- dplyr::filter(s_QTL, Chr == p)
      
      # Sort the regions based on their start positions
      temp <- temp[order(temp$Start), ]
      
      # Initialize variables to keep track of the best region convergence
      best_start <- 0
      best_end <- 0
      max_convergence <- 1
      
      # Conditional to see if there is only one QTL per chromosome 
      if (nrow(temp) > 1){
        
        # Iterate over each region and find the convergence
        for (i in 1:(nrow(temp) - 1)){
          
          # Get the current region
          region_start <- temp$Start[i]
          region_end <- temp$End[i]
          
          # Calculate the convergence of the current region
          convergence <- 1
          
          # Iterate over the subsequent regions
          for (j in (i + 1):nrow(temp)){
            
            # Get the next region
            next_region_start <- temp$Start[j]
            next_region_end <- temp$End[j]
            
            # Check if the next region is within the current region
            if (next_region_start >= region_start && next_region_end <= region_end){
              
              # Update the convergence count
              convergence <- convergence + 1
              
              # Update region start
              region_start <- next_region_start
              
            } else {
              
              # Save the best convergence regions
              best_start <- region_start
              best_end <- region_end
              
            }
          }
          
          # Check if the current region has higher convergence than the previous best
          if (convergence > max_convergence){
            best_start <- region_start
            best_end <- region_end
            max_convergence <- convergence
            
          }
        }
        
      } else {
        
        # Gives the values directly
        best_start <- temp$Start[1]
        best_end <- temp$End[1]
        Convergence <- 1
        Chr <- p
        
      }
      
      # Print the region with the highest convergence
      best[[p]] <- data.frame(Chr = p, Start = best_start, End = best_end, Convergence = max_convergence)
      
    }
    
    # Count the number of QTLs per chromosome
    num_chrs <- data.frame(table(s_QTL$Chr))
    colnames(num_chrs) <- c("Chr", "QTL_count")
    num_chrs$Chr <- as.integer(as.character(num_chrs$Chr))
    
    # Merge the data
    c_QTL <- bind_rows(best)
    c_QTL <- inner_join(c_QTL, num_chrs, by = "Chr")
    c_QTL <- c_QTL %>% mutate(QTl.Wide = End - Start) %>%
      rename(QTL.Start = Start, QTL.End = End, QTL.Convergence = Convergence) %>%
      select(Chr, QTL.Start, QTL.End, QTL.Convergence, QTL_count, QTl.Wide)
    
    # Informative message
    message("QTL files loaded succesfully")
    
    } else {
    
    # Set working directory
    setwd(Wdir)
    
    # Read the files with the results
    # Creates the single QTL data frame (s_QTL)
    s_QTL <- read.csv(paste0(name)) %>%
      select(phenotype, chr, lod, start.marker, end.marker) %>%
      mutate(start = start.marker) %>%
      tidyr::separate(col = start, into = c("na", "Start"), sep = paste("_"), extra = "drop") %>%
      mutate(end = end.marker) %>%
      tidyr::separate(col = end, into = c("na", "End"), sep = paste("_"), extra = "drop") %>%
      select(phenotype, chr, lod, start.marker, end.marker, Start, End) %>%
      rename(Phenotype = phenotype, Chr = chr, LOD = lod,
             Start.Marker = start.marker, End.Marker = end.marker)
    
    # Re format columns
    s_QTL$Start <- as.numeric(s_QTL$Start)
    s_QTL$End <- as.numeric(s_QTL$End)
    
    
    
    # Creates the converged QTL data frame (c_QTL)
    # Creates am empty list to save the QTLs converge regions per chromosomes 
    best <- list()
    
    # Find the regions in which converge most of QTLs per chromosome
    for (p in unique(s_QTL$Chr)){
      
      # Filter the QTL data per chromosome
      temp <- dplyr::filter(s_QTL, Chr == p)
      
      # Sort the regions based on their start positions
      temp <- temp[order(temp$Start), ]
      
      # Initialize variables to keep track of the best region convergence
      best_start <- 0
      best_end <- 0
      max_convergence <- 1
      
      # Conditional to see if there is only one QTL per chromosome 
      if (nrow(temp) > 1){
        
        # Iterate over each region and find the convergence
        for (i in 1:(nrow(temp) - 1)){
          
          # Get the current region
          region_start <- temp$Start[i]
          region_end <- temp$End[i]
          
          # Calculate the convergence of the current region
          convergence <- 1
          
          # Iterate over the subsequent regions
          for (j in (i + 1):nrow(temp)){
            
            # Get the next region
            next_region_start <- temp$Start[j]
            next_region_end <- temp$End[j]
            
            # Check if the next region is within the current region
            if (next_region_start >= region_start && next_region_end <= region_end){
              
              # Update the convergence count
              convergence <- convergence + 1
              
              # Update region start
              region_start <- next_region_start
              
            } else {
              
              # Save the best convergence regions
              best_start <- region_start
              best_end <- region_end
              
            }
          }
          
          # Check if the current region has higher convergence than the previous best
          if (convergence > max_convergence){
            
            # Save the best convergence regions
            best_start <- region_start
            best_end <- region_end
            max_convergence <- convergence
            
          }
        }
        
      } else {
        
        # Gives the values directly
        best_start <- temp$Start[1]
        best_end <- temp$End[1]
        Convergence <- 1
        Chr <- p
        
      }
      
      # Print the region with the highest convergence
      best[[p]] <- data.frame(Chr = p, Start = best_start, End = best_end, Convergence = max_convergence)
      
    }
    
    # Count the number of QTLs per chromosome
    num_chrs <- data.frame(table(s_QTL$Chr))
    colnames(num_chrs) <- c("Chr", "QTL_count")
    num_chrs$Chr <- as.integer(as.character(num_chrs$Chr))
    
    # Merge the data
    c_QTL <- bind_rows(best)
    c_QTL <- inner_join(c_QTL, num_chrs, by = "Chr")
    c_QTL <- c_QTL %>% mutate(QTl.Wide = End - Start) %>%
      rename(QTL.Start = Start, QTL.End = End, QTL.Convergence = Convergence) %>%
      select(Chr, QTL.Start, QTL.End, QTL.Convergence, QTL_count, QTl.Wide)
    
    # Informative message
    message("QTL files loaded succesfully")
    
    }
  
  
  
  # 3: Match the QTL results with the gene annotation databases ----------------
  
  # Conditional: If there is at least one QTL to annotate, the function continues
  if (dim(s_QTL)[1] > 0){
    
    # Informative message
    message(paste("There are", dim(s_QTL)[1], "QTLs to annotate"))
    
    
    
    # 3.1 Function pathway for single QTLs (s_QTL) -----------------------------
    
    # 3.1.1 Starts filtering the gff3 data -------------------------------------
    # Informative message
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
      cat('\r', i, ' rows processed |', rep('=', i / 4), ifelse(i == nrow(s_QTL), '|\n', '>'), sep = '')
      
    }
    
    # Merge all the data frames of the list in a unique data frame
    s_gff_m <- bind_rows(s_gff_l) %>%
      select(Trait, Chr, LOD, QTL.Start, QTL.End, What, Name, Location, Gen.Start, Gen.End)
    
    
    
    # 3.2.2 Continues filtering the annotation data ----------------------------
    # Informative message
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
        mutate(Trait = s_gff_m$Trait[i], Chr = s_gff_m$Chr[i], LOD = s_gff_m$LOD[i],
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
    
    
    
    # 3.2 Function pathway for merged QTLs (c_QTL) -----------------------------
    
    # 3.2.1 Starts filtering the gff3 data -------------------------------------
    # Informative messages
    message("Obtaining gene information...")
    
    # Creates an empty list object to save the data
    c_gff_l <- list()
    
    # Filters and merges the gff3 data
    for (i in 1:nrow(c_QTL)){
      
      # Filter the gff3 file per chromosome
      data <- dplyr::filter(gff, c_QTL$Chr[i] == Chr)
      
      # Combine the different filtered data
      c_gff_l[[i]] <- rbind(
        
        # Looking for genes located inside the QTL
        dplyr::filter(data, c_QTL$QTL.Start[i] <= Start & c_QTL$QTL.Start[i] <= End &
                        c_QTL$QTL.End[i] >= Start & c_QTL$QTL.End[i] >= End) %>%
          mutate(QTL.Start = c_QTL$QTL.Start[i], QTL.End = c_QTL$QTL.End[i],
                 Location = "Gene inside QTL") %>%
          rename(Gen.Start = Start, Gen.End = End),
        
        # Looking for genes that contain the QTL
        dplyr::filter(data, c_QTL$QTL.Start[i] >= Start & c_QTL$QTL.Start[i] <= End &
                        c_QTL$QTL.End[i] >= Start & c_QTL$QTL.End[i] <= End) %>%
          mutate(QTL.Start = c_QTL$QTL.Start[i], QTL.End = c_QTL$QTL.End[i],
                 Location = "QTL inside gene") %>%
          rename(Gen.Start = Start, Gen.End = End),
        
        # Looking genes down stream 
        dplyr::filter(data, c_QTL$QTL.Start[i] >= Start & c_QTL$QTL.Start[i] <= End) %>%
          mutate(QTL.Start = c_QTL$QTL.Start[i], QTL.End = c_QTL$QTL.End[i],
                 Location = "Gene overlaps QTL") %>%
          rename(Gen.Start = Start, Gen.End = End),
        
        # Looking genes up stream
        dplyr::filter(data, c_QTL$QTL.End[i] >= Start & c_QTL$QTL.End[i] <= End) %>%
          mutate(QTL.Start = c_QTL$QTL.Start[i], QTL.End = c_QTL$QTL.End[i],
                 Location = "Gene overlaps QTL") %>%
          rename(Gen.Start = Start, Gen.End = End)
        
      )
      
      # Progress bar
      cat('\r', i, ' files processed |', rep('=', i / 1), ifelse(i == nrow(c_QTL), '|\n', '>'), sep = '')
      
    }
    
    # Merge the data frames of the list in a single data frame
    c_gff_m <- bind_rows(c_gff_l) %>%
      select(Chr, QTL.Start, QTL.End, What, Name, Location, Gen.Start, Gen.End)
    
    
    
    # 3.2.2 Continues filtering the annotation data ----------------------------
    # Informative message
    message("Obtaining gene annotation...")
    
    # Creates an empty list object to save the data
    c_annot_l <- list()
    
    # Filters and merges the annot data
    for (i in 1:nrow(c_gff_m)){
      
      # Filter the annotation data and creates new columns based on previous data
      c_annot_l[[i]] <- filter(annot,
                                Locus == c_gff_m$Name[i] |
                                  Trans == c_gff_m$Name[i] |
                                    Peptide == c_gff_m$Name[i]) %>%
        mutate(Chr = c_gff_m$Chr[i], What = c_gff_m$What[i], Location = c_gff_m$Location[i],
               Gen.Start = c_gff_m$Gen.Start[i], Gen.End = c_gff_m$Gen.End[i],
               QTL.Start = c_gff_m$QTL.Start[i], QTL.End = c_gff_m$QTL.End[i])
      
      # Progress bar
      cat('\r', i, ' rows processed |', rep('=', i / 10), ifelse(i == nrow(c_annot_l), '|\n', '>'), sep = '')
      
    }
    
    
    
    # 3.2.3 It ends by formatting the data frame -------------------------------
    
    # Merge the data frames of the list in a single data frame and modify it
    c_QTL_annotation <- bind_rows(c_annot_l) %>% 
      select(Chr, QTL.Start, QTL.End, What, ID, Locus, Trans, Peptide, Location, 
             Gen.Start, Gen.End, GO, AT.name, AT.define)
    
    
    
    # 4: Saving the outputs ----------------------------------------------------
    
    # Create empty lists
    s_QTL_annotation_l <- list()
    c_QTL_annotation_l <- list()
    
    # Split the data in a list by saving each chromosome in a different data frame
    for (p in sort(unique(s_QTL$Chr))){
      
      # Filter the QTL data per chromosome
      s_QTL_annotation_l[[p]] <- dplyr::filter(s_QTL_annotation, Chr == p)
      c_QTL_annotation_l[[p]] <- dplyr::filter(c_QTL_annotation, Chr == p)
      
    }
    
    # Create empty Excel workbooks
    s_QTL_annotation_wb <- createWorkbook()
    c_QTL_annotation_wb <- createWorkbook()
    
    # Add sheets to the workbook
    for (p in sort(unique(s_QTL$Chr))){
      
      # Set sheet name
      addWorksheet(s_QTL_annotation_wb, sheetName = paste0("Chr_0", p))
      addWorksheet(c_QTL_annotation_wb, sheetName = paste0("Chr_0", p))
      
      # Obtain data
      writeData(s_QTL_annotation_wb, sheet = paste0("Chr_0", p),
                x = s_QTL_annotation_l[[p]], startCol = 1, startRow = 1)
      writeData(c_QTL_annotation_wb, sheet = paste0("Chr_0", p),
                x = c_QTL_annotation_l[[p]], startCol = 1, startRow = 1)
      
    }
    
    # Save the workbooks
    # Informative messages
    message("Saving output files as: 'single_QTL_annotation.xlsx' and 'merged_QTL_annotation.xlsx'")
    
    saveWorkbook(s_QTL_annotation_wb, paste0(prefix, ".single_QTL_annotation.xlsx"), overwrite = T)
    saveWorkbook(c_QTL_annotation_wb, paste0(prefix, ".merged_QTL_annotation.xlsx"), overwrite = T)
    
    message("Done!")
    
  } else {
    
    # Informative message
    message("No QTLs to annotate")
    
  }
  
}



# Example(s) -------------------------------------------------------------------
# Set arguments
# Ddir <- "D:/OneDrive - CGIAR/00_CassavaBioinformaticsPlatform/00_Basics/00_Data/Mesculenta_305_v6.1/"
# annot <- "Mesculenta_305_v6.1.annotation_info.txt"
# gff <- "Mesculenta_305_v6.1.gene.gff3"
# version <- "6.1"
# wdyw <- "gene"
# prefix <- "Candiate_genes_WFR"

# If recursive
# Wdir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F1_Metabolomics/02_QTL_Analysis/CM8996/Everything/"
# name <- "LodIntervals"
# recursive <- "T"

# Non-recursive
# Wdir <- "D:/OneDrive - CGIAR/00_CassavaBioinformaticsPlatform/01_ACWP/08_CandidateGenes/"
# name <- "RegionQTL.csv"
# recursive <- "F"



# Run function -----------------------------------------------------------------
# QTL_Annotation(Ddir, annot, gff, version, wdyw, prefix, Wdir, name, recursive)
