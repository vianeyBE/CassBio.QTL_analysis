# Short name: Annotation for QTL results
# Description (1): Get the annotation of the gen containing the SNP using the rqtl/QTL maping results
# Description (2): Retrieves the closest genes down and up stream
# Output: Basename of the output file
#
# Authors: Camilo E. Sanchez (c.e.sanchez@cgiar.org) and Vianey Barrera-Enriquez (vpbarrera@gmail.com)
#
# Arguments:
# Mdir: Name of the directory that contains the GAPIT results. For example: home/user/folder.
# pat: Enter the path of file names to look for. For example: QTL_LOD_Intervals. The path must finish with a point (.).
# wdyw: Enter what are you looking for to annotate. Options: CDS, five_prime_UTR, gene, mRNA, three_prime_UTR.



###### To do ######
# 1. Transform the script into a function



###### Examples ######
Mdir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/03_GWAS_PPD_Populations/04_GWAS/GAPIT_Results"
pat <- "GAPIT.Association.GWAS_Results."
wdyw <- "gene"



# 0: Function init -------------------------------------------------------------

QTL_Annotation <- function(Mdir, pat, wdyw){
  
  # 1: Load all the directories, info, and data --------------------------------
  
  # Load packages
  if (!require(tidyverse)) install.packages(tidyverse)
  library(tidyverse)
  
}

setwd("D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/00_Data/Mesculenta_305_v6.1")

message("Loading required files to do the annotation")

annot <- read.delim("Mesculenta_305_v6.1.annotation_info.txt", header = F) %>%
  rename(ID = V1, Gen1 = V2, Gen2 = V3, Gen3 = V4, PF = V5, PTHR = V6, KOG = V7, na = V8, K0 = V9,
         GO = V10, AT = V11, NPI = V12, name = V13) %>%
  select(ID, Gen1, Gen2, Gen3, GO, NPI, name)

gff3 <- read.delim("Mesculenta_305_v6.1.gene.gff3", header = F, comment.char = "#") %>%
  rename(chr = V1, phyto = V2, what = V3, start = V4, end = V5, na = V6, sign = V7, NPI = V8, sep = V9) %>%
  tidyr::separate(col = sep, into = c("ID", "na"), sep = ";") %>%
  separate(col = ID, into = c("na2", "na3"), sep = "=") %>%
  separate(col = na3, into = c("name", "na4"), sep = ".v" ) %>%
  select(chr, what, start, end, name) %>%
  mutate(chr = recode(chr, Chromosome01 = 1, Chromosome02 = 2, Chromosome03 = 3, Chromosome04 = 4,
                      Chromosome05 = 5, Chromosome06 =  6, Chromosome07 = 7, Chromosome08 = 8, 
                      Chromosome09 = 9, Chromosome10 = 10, Chromosome11 = 11, Chromosome12 = 12,
                      Chromosome13 = 13, Chromosome14 = 14, Chromosome15 = 15, Chromosome16 = 16, 
                      Chromosome17 = 17, Chromosome18 = 18)) %>%
  na.omit() %>% 
  dplyr::filter(what %in% wdyw)



# 2: Find all the CSVs with the results and filter them ------------------------

message("Getting list of CSV files...")

# Get the names of the files
setwd(Mdir)
names <- list.files(path = Mdir, pattern = pat, all.files = F, full.names = F, recursive = T)

message(paste("QTL mapping files found:", length(names)))

# Create an empty list and the variable to iterate
# Read, rename, modify and save all the csv files in the list 
message("Reading QTL mapping files...")

# Creates the objects to perform the loop
csvL <- list()
p <- 0

# Loop
for (i in 1:length(names)){
  p <- p + 1
  
  name.F <- gsub(" / ", "/", paste(Mdir, "/", names[p]))
  name.T <- gsub(paste(Mdir), "", paste(name.F))
  
  csvL[[i]] <- read.delim(name.F) %>% select(phenotype, chr, lod, start.marker, end.marker) %>%
    mutate(start = start.marker) %>%
    tidyr::separate(col = start, into = c("na", "start"), sep = paste("_")) %>%
    mutate(end = end.marker) %>%
    tidyr::separate(col = end, into = c("na", "end"), sep = paste("_")) %>%
    select(phenotype, chr, lod, start.marker, end.marker, start, end)
  
  # Progress bar
  cat('\r', i, ' files processed |', rep('=', i / 4), ifelse(i == length(names), '|\n', '>'), sep = '')
}

# Merge the data frames of the list in a single data frame
QTL <- bind_rows(csvL)
QTL$start <- as.numeric(QTL$start)
QTL$end <- as.numeric(QTL$end)
QTL_s <- QTL %>% group_by(chr) %>% summarise("QTL_start" = min(start), "QTL_end" = max(end))

message(paste("There are", dim(QTL)[1], "QTLs after filtering"))

# Delete unnecessary objects
rm(csvL, i, p, name.F, name.T)



# 3: Match the results with the gene annotation database -----------------------

# Creates a conditional
# If there is at least one result of the GWAS models, the annotation is done
if (dim(QTL_s)[1] > 0){
  
  message("Retriving annotation...")
  
  # Creates the objects to perform the loop
  p <- 0
  annotL <- list()
  
  # Loop
  for (i in 1:nrow(QTL_s)){
    p <- p + 1
    
    data <- dplyr::filter(gff3, QTL_s$chr[p] == chr)
    
    annotL[[i]] <- rbind(
      # Looking for genes located inside the QTL
      dplyr::filter(data, QTL_s$QTL_start[p] <= start & QTL_s$QTL_start[p] <= end &
                    QTL_s$QTL_end[p] >= start & QTL_s$QTL_end[p] >= end) %>%
      mutate(QTL_start = QTL_s$QTL_start[p], QTL_end = QTL_s$QTL_end[p], Location = "Gene inside QTL") %>%
      rename(gen_start = start, gen_end = end),
      
      # Looking for genes that contain the QTL
      dplyr::filter(data, QTL_s$QTL_start[p] >= start & QTL_s$QTL_start[p] <= end &
                    QTL_s$QTL_end[p] >= start & QTL_s$QTL_end[p] <= end) %>%
      mutate(QTL_start = QTL_s$QTL_start[p], QTL_end = QTL_s$QTL_end[p], Location = "QTL inside gene") %>%
      rename(gen_start = start, gen_end = end),
      
      # Looking genes down stream 
      dplyr::filter(data, QTL_s$QTL_start[p] >= start & QTL_s$QTL_start[p] <= end) %>%
      mutate(QTL_start = QTL_s$QTL_start[p], QTL_end = QTL_s$QTL_end[p], Location = "Gene overlaps QTL") %>%
      rename(gen_start = start, gen_end = end),
      
      # Looking genes up stream
      dplyr::filter(data, QTL_s$QTL_end[p] >= start & QTL_s$QTL_end[p] <= end) %>%
      mutate(QTL_start = QTL_s$QTL_start[p], QTL_end = QTL_s$QTL_end[p], Location = "Gene overlaps QTL") %>%
      rename(gen_start = start, gen_end = end)
      )
    
    # Progress bar
    cat('\r', i, ' files processed |', rep('=', i / 4), ifelse(i == nrow(QTL_s), '|\n', '>'), sep = '')
    }

  # Merge the data frames of the list in a single data frame
  annotLm <- bind_rows(annotL)

  message(paste("There are", dim(annotLm)[1], "genes that overlap with QTLs"))

  # Delete unnecessary objects
  rm(annotL, data, i, p)



  # 4: Formatting dataframe ----------------------------------------------------

  message("Obtaining gene information...")
  
  # Creates the objects to perform the loop
  p <- 0
  gffL <- list()
    
  # Loop
  for (i in 1:nrow(annotLm)){
    p <- p + 1
    
    gffL[[i]] <- filter(annot, Gen1 == annotLm$name[p] | Gen2 == annotLm$name[p] | Gen3 == annotLm$name[p]) %>%
       mutate(chr = annotLm$chr[p], 
              gen_start = annotLm$gen_start[p], gen_end = annotLm$gen_end[p],
              QTL_start = annotLm$QTL_start[p], QTL_end = annotLm$QTL_end[p], 
              Location = annotLm$Location[p])
      
    # Progress bar
    cat('\r', i, ' files processed |', rep('=', i / 4), ifelse(i == nrow(annotLm), '|\n', '>'), sep = '')
    }

    # Merge the data frames of the list in a single data frame and modify it
    QTL_annotation <- bind_rows(gffL) %>% 
      select(chr, Gen1, name, gen_start, gen_end, GO, NPI, QTL_start, QTL_end, Location) %>%
      rename(gen_name = Gen1, gen_name_extend = name)
    
    rm(annotLm, gffL, i, p)
    
    # 5: Save output -----------------------------------------------------------
    
    message("Saving output file: 'QTL_Annotation.csv' in working directory")
    
    write.csv(QTL_annotation, file = "QTL_Annotation.csv", quote = F, row.names = F)
    
    message("Done! ")
    
} else {
  
  message("No QTLs to annotate")
  
  }
