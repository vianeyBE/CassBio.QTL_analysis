# Short name: Annotation for QTL results
# Description (1): Get the annotation of region using the QTL results
# Description (2): Retrieves genes that overlap or are inside the region
# Output: Excel files with the annotation results
# 
# Authors (1): Camilo E. Sanchez (c.e.sanchez@cgiar.org)
# Authors (2): Vianey Barrera-Enriquez (vpbarrera@gmail.com)
# 
# Arguments:
# dir: Path to the folder that contains the QTL results (e.g., /home/user/folder)
# name: Name of file names to look for (e.g., QTL_Intervals)
# wdyw: What to annotate (Opts: CDS, five_prime_UTR, gene, mRNA, three_prime_UTR)
# recursive: Perform a recursive search or not (Boolean string. Default = F)
# version: Genome of reference version (Opts: 6.1 or 8.1)
# prefix: Prefix name. How the outputs will be saved (e.g., Flowering)




# To do ------------------------------------------------------------------------
# 1: Add trait name in c_QTL 




# Example ----------------------------------------------------------------------
# Set arguments
# dir <- "D:/OneDrive - CGIAR/00_BioInf_Platform/16_Flowering/02_Annotation/"
# name <- "Regions_QTLs.csv"
# recursive <- FALSE
# version <- "6.1"
# wdyw <- "gene"
# prefix <- "Flowering"




# 0: Function init -------------------------------------------------------------
QTL_Annotation <- function(dir, name, recursive = F, version, wdyw, prefix) {
  
  
  
  
  # 1: Set working directory and load package ----------------------------------
  # Load packages
  if (!require("tidyverse")) install.packages("tidyverse")
  if (!require("openxlsx")) install.packages("openxlsx")
  
  library(tidyverse)
  library(openxlsx)
  
  
  
  
  # 2: Build helper functions --------------------------------------------------
  # Helper 1: Load annotation data (annot and gff3 files)
  load_annotation <- function(version = c("6.1", "8.1"), wdyw) {
    
    # Verify arguments
    version <- match.arg(as.character(version), c("6.1", "8.1"))
    
    # Informative message
    message("Loading required files to do the annotation...")
    
    # Set the path for the data
    dir_61 <- "D:/OneDrive - CGIAR/00_BioInf_Platform/00_Basics/03_Server_Bioinfo-CO/00_ReferenceGenome/Mesculenta_v6.1/"
    dir_81 <- "D:/OneDrive - CGIAR/00_BioInf_Platform/00_Basics/03_Server_Bioinfo-CO/00_ReferenceGenome/Mesculenta_v8.1/"
    
    # Resolve paths by version 
    if (version == "6.1") {
      root        <- dir_61
      annot_file  <- "Mesculenta_305_v6.1.annotation_info.txt"
      gff_file    <- "Mesculenta_305_v6.1.gene.gff3"
      has_header  <- FALSE
    } else {
      root        <- dir_81
      annot_file  <- "Mesculenta_671_v8.1.annotation_info.txt"
      gff_file    <- "Mesculenta_671_v8.1.gene.gff3"
      has_header  <- TRUE
    }
    
    # Set data paths
    annot_path <- file.path(root, annot_file)
    gff_path   <- file.path(root, gff_file)
    
    # Check if the files exist
    if (!file.exists(annot_path)) stop("Annotation file not found: ", annot_path)
    if (!file.exists(gff_path))   stop("GFF file not found: ", gff_path)
    
    # Read annotation table
    if (version == "6.1") { # v6.1: no header; keep ID, Name, GO, Annotation
      annot_raw <- readr::read_tsv(annot_path, col_names = F, show_col_types = F)
      annot <- tibble(Name       = annot_raw[[2]],
                      GO         = annot_raw[[10]],
                      Annotation = annot_raw[[13]])
    } else {
      # v8.1: header present; keep Name, GO, Annotation (ID not provided)
      annot_raw <- readr::read_tsv(annot_path, col_names = T, show_col_types = F)
      annot <- tibble(Name       = annot_raw[[2]],
                      GO         = annot_raw[[10]],
                      Annotation = annot_raw[[12]])
    }
    
    # Read and tidy GFF3
    # Columns: Chr (V1), What (V3), Start (V4), End (V5), Strand (V7), Attr (V9)
    gff_raw <- readr::read_tsv(gff_path, col_names = F, comment = "#", show_col_types = F)
    
    gff <- gff_raw %>% transmute(Chr       = X1,
                                 What      = X3,
                                 Gen.Start = as.numeric(X4),
                                 Gen.End   = as.numeric(X5),
                                 Strand    = X7,
                                 attr      = X9) %>%
      mutate(Name = stringr::str_match(attr, "ID=([^;]+)")[, 2],
             Name = stringr::str_remove(Name, "\\.v.*$")) %>%
      select(Chr, What, Gen.Start, Gen.End, Name, Strand) %>%
      filter(stringr::str_detect(Chr, "Chromosome")) %>%
      mutate(Chr = as.integer(stringr::str_extract(Chr, "\\d+"))) %>%
      filter(What %in% wdyw)
    
    # Informative message
    message("Annotation files loaded successfully")
    
    # Retrieve data
    list(annot = annot, gff = gff)
    
  }
  
  # Helper 2: Read and standardize a single QTL file
  parse_qtl <- function(path) {
    
    readr::read_csv(path, show_col_types = FALSE) %>%
      select(any_of(c("phenotype", "chr", "lod", "start.marker", "end.marker"))) %>%
      rename(Phenotype    = phenotype, 
             Chr          = chr,
             Start.Marker = start.marker,
             End.Marker   = end.marker) %>% # Standardize names (keep 'lod' if present)
      mutate(QTL.Start = as.numeric(stringr::str_extract(Start.Marker, "(?<=_)\\d+(?:\\.\\d+)?")),
             QTL.End   = as.numeric(stringr::str_extract(End.Marker, "(?<=_)\\d+(?:\\.\\d+)?")),
             Chr       = as.character(Chr)) %>% # Extract numeric position after the last underscore
      select(Phenotype, Chr, Start.Marker, End.Marker, QTL.Start, QTL.End, any_of("lod"))
    
  }
  
  # Helper 3: Find per-chromosome convergence with a sweep-line algorithm
  converge_qtl <- function(s_QTL) {
    
    if (nrow(s_QTL) == 0) {
      return(tibble(Chr             = character(),
                    QTL.Start       = numeric(),
                    QTL.End         = numeric(),
                    QTL.Convergence = integer(),
                    QTL_count       = integer(),
                    QTl.Wide        = numeric()))
    
    }
    
    # Per-chromosome: segment with maximum overlap (ties -> widest, then earliest start)
    c_QTL <- s_QTL %>% group_by(Chr) %>%
      group_modify(~{
        
        # 
        df <- .x
        n  <- nrow(df)
        
        # 
        if (n == 1) {
          
          #
          tibble(QTL.Start = df$QTL.Start, QTL.End = df$QTL.End, QTL.Convergence = 1L)
          
        } else {
          
          #
          ev <- bind_rows(tibble(pos = df$QTL.Start, delta =  1L),
                          tibble(pos = df$QTL.End,   delta = -1L)) %>%
            arrange(pos, desc(delta)) # starts before ends at the same position
          
          #
          k <- cumsum(ev$delta)
          
          #
          segs <- tibble(start = ev$pos[-length(ev$pos)],
                         end = ev$pos[-1],
                         cov = k[-length(k)]) %>%
            filter(end > start)
          
          #
          if (nrow(segs) == 0) {
            #
            tibble(QTL.Start       = min(df$QTL.Start),
                   QTL.End         = max(df$QTL.End),
                   QTL.Convergence = n)
          } else {
            
            #
            max_cov <- max(segs$cov)
            
            # 
            segs %>%
              filter(cov == max_cov) %>%
              mutate(width = end - start) %>%
              slice_max(width, n = 1, with_ties = TRUE) %>%
              slice_min(start, n = 1, with_ties = TRUE) %>%
              transmute(QTL.Start       = start, 
                        QTL.End         = end, 
                        QTL.Convergence = as.integer(cov))
            
          }
          
        }
        
      }) %>%
      ungroup() %>%
      left_join(count(s_QTL, Chr, name = "QTL_count"), by = "Chr") %>%
      mutate(QTl.Wide = QTL.End - QTL.Start) %>%
      select(Chr, QTL.Start, QTL.End, QTL.Convergence, QTL_count, QTl.Wide)
    
    # Retrieve data
    c_QTL
    
  }
  
  # Helper 4: Unified loader 
  load_qtl <- function(dir, name, recursive = FALSE) {
    
    # 
    if (isTRUE(recursive)) {
      
      # Informative message
      message("Getting list of CSV files...")
      
      # Pattern: files that start with 'name' and end with .csv (adjust if needed)
      files <- list.files(path = dir,
                          pattern = paste0("^", name, ".*\\.csv$"),
                          full.names = TRUE,
                          recursive = TRUE)
      
      # Informative message
      message(paste("QTL mapping files found:", length(files)))
      
      #
      if (length(files) == 0) stop("No files found. Check 'dir', 'name' and extension.")
      
      # Informative message
      message("Reading QTL mapping files...") 
      
    } else {
      
      # Informative message
      message("Loading QTL file...")
      
      #
      files <- name
      
      #
      if (!file.exists(files)) stop("File not found: ", files)
      
    }
    
    # s_QTL: single long table of all QTLs
    s_QTL <- purrr::map_dfr(files, parse_qtl) %>%
      mutate(Chr = as.numeric(Chr),
             QTL.Start = as.numeric(QTL.Start),
             QTL.End = as.numeric(QTL.End))
    
    # c_QTL: per-chromosome convergence
    c_QTL <- converge_qtl(s_QTL) %>%
      mutate(Chr = as.numeric(Chr),
             QTL.Start = as.numeric(QTL.Start),
             QTL.End = as.numeric(QTL.End))
    
    # Informative message
    message("QTL file(s) loaded successfully")
    
    # Retrieve data
    list(s_QTL = s_QTL, c_QTL = c_QTL)
    
  }
  
  # Helper 5: Classify gene–QTL relationship
  classify_location <- function(df) {
    df %>% mutate(Location = case_when(
      Gen.Start >= QTL.Start & Gen.End <= QTL.End ~ "Gene inside QTL",
      QTL.Start >= Gen.Start & QTL.End <= Gen.End ~ "QTL inside gene",
      TRUE                                        ~ "Gene overlaps QTL"))
  }
  
  # Helper 6: # Filter gene info and match it with QTL data
  annotate_genes <- function(s_QTL, c_QTL, gff) {
    
    # Informative message
    message("Matching genes to s_QTL and c_QTL intervals...")
    
    # s_QTL x GFF
    s_gff_m <- gff %>%
      inner_join(s_QTL %>% transmute(Chr, QTL.Start, QTL.End, Trait = Phenotype),
                 by = "Chr") %>%
      filter(Gen.Start <= QTL.End, Gen.End >= QTL.Start) %>%
      classify_location() %>%
      select(Trait, Chr, QTL.Start, QTL.End, What, Name, Location, Gen.Start, Gen.End, Strand) %>%
      arrange(Chr, QTL.Start, Gen.Start)
    
    # c_QTL x GFF
    c_gff_m <- gff %>%
      inner_join(c_QTL %>% transmute(Chr, QTL.Start, QTL.End),
                 by = "Chr") %>%
      filter(Gen.Start <= QTL.End, Gen.End >= QTL.Start) %>%
      classify_location() %>%
      select(Chr, QTL.Start, QTL.End, What, Name, Location, Gen.Start, Gen.End, Strand) %>%
      arrange(Chr, QTL.Start, Gen.Start)
    
    # Informative message
    message(paste("There are", nrow(s_gff_m), "genes to annotate (s_QTL)"))
    message(paste("There are", nrow(c_gff_m), "genes to annotate (c_QTL)"))
    
    # Return both annotated tables
    list(s_gff = s_gff_m, c_gff = c_gff_m)
    
  }
  
  # # Helper 7: build a named list of data frames per chromosome
  sheet_list <- function(df, chr_vec, label = "Chr") {
    
    # Empty frame with same columns
    proto <- df[0, ]
    
    # 
    setNames(lapply(chr_vec, function(ch) {
      
      # 
      tmp <- dplyr::filter(df, Chr == ch)
      
      #
      if (nrow(tmp) == 0) proto else tmp}),
      
      sprintf("%s_%02d", label, as.integer(chr_vec)))
    
  }
  
  
  
  
  # 3: Execute the helpers to find and load the data ---------------------------
  # Set working directory
  setwd(dir)
  
  # Execute function to load annotations files)
  annotations <- load_annotation(version, wdyw)
  
  # Execute function to load QTL file(s)
  qtls <- load_qtl(dir, name, recursive)
  
  # Extract the dataframes
  annot <- annotations$annot
  gff   <- annotations$gff
  s_QTL <- qtls$s_QTL
  c_QTL <- qtls$c_QTL
  
  # Informative message
  message(paste0("There are ", dim(s_QTL)[1], " single QTLs and ", 
                 dim(s_QTL)[1], " convergence QTLs to annotate"))
  
  
  
  
  # 4: Filtering the gff3 data -------------------------------------------------
  gff_clean <- annotate_genes(s_QTL, c_QTL, gff)
  s_gff <- gff_clean$s_gff
  c_gff <- gff_clean$c_gff
    
    
  
    
  # 5: Match the annotation data -----------------------------------------------
  # Informative messages
  message("Joining annotation for single QTLs (s_QTL) and",
          "convergence QTLs (c_QTL) genes...")
  
  # Vectorized join for s_QTl
  s_QTL_annotation <- s_gff %>%
    inner_join(annot, by = "Name") %>%
    mutate(ID = paste0(Chr, Name, Gen.Start, Gen.End, Strand)) %>%
    filter(!duplicated(ID)) %>%
    select(-ID, Trait, Chr, QTL.Start, QTL.End, What, Name, Location, 
           Gen.Start, Gen.End, Strand, GO, Annotation) %>%
    arrange(Chr, QTL.Start, Gen.Start, Name)
  
  # Vectorized join for c_QTL (keeps only genes with annotation)
  c_QTL_annotation <- c_gff %>%
    inner_join(annot, by = "Name") %>%
    mutate(ID = paste0(Chr, Name, Gen.Start, Gen.End, Strand)) %>%
    filter(!duplicated(ID)) %>%
    select(-ID, Chr, QTL.Start, QTL.End, What, Name, Location, 
           Gen.Start, Gen.End, Strand, GO, Annotation) %>%
    arrange(Chr, QTL.Start, Gen.Start, Name)
  
  # Informative messages
  message(paste0("There are", nrow(s_QTL_annotation), 
                 " genes annotated for single QTLs (s_QTL) and ",
                 nrow(c_QTL_annotation), 
                 " genes annotated for convergence QTLs (c_QTL)"))
    
    
    
    
  # 6: Saving outputs -----------------------------------------------------------
  message("Preparing Excel workbooks...")
  
  # Decide sheet order
  chr_order <- if (exists("s_QTL", inherits = FALSE)) {
    sort(unique(s_QTL$Chr))
  } else {
    sort(unique(c(s_QTL_annotation$Chr, c_QTL_annotation$Chr)))
  }
  
  # Build sheet lists (one per workbook)
  s_sheets <- sheet_list(s_QTL_annotation, chr_order)
  c_sheets <- sheet_list(c_QTL_annotation, chr_order)
  
  # Informative message
  message(paste0("\nSaving output files as:\n",
                 "- ", prefix, "_sQTL.annotation.xlsx\n",
                 "- ", prefix, "_cQTL.annotation.xlsx"))
  
  # Write both Excel files (one sheet per chromosome)
  openxlsx::write.xlsx(x         = s_sheets,
                       file      = paste0(prefix, "_sQTL.annotation.xlsx"),
                       overwrite = TRUE)
  
  openxlsx::write.xlsx(x         = c_sheets,
                       file      = paste0(prefix, "_cQTL.annotation.xlsx"),
                       overwrite = TRUE)
  
  # Informative message
  message("All is done!")
  
}




# Run function -----------------------------------------------------------------
# QTL_Annotation(dir, name, recursive, version, wdyw, prefix)
