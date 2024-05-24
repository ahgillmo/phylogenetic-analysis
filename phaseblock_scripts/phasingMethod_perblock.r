#!/usr/bin
# Required packages
suppressWarnings(library(dplyr, warn.conflicts = FALSE))

# Input command is the path to the input file in this case it would be the bc.txt files
input=commandArgs(trailing=TRUE)

# Read in BC files
BCfile <- read.table(file = input[1], sep = "", header = FALSE, stringsAsFactors = FALSE)

# Remove the duplicated variants
BCfile <- BCfile[!duplicated(BCfile$V1), ]

# Phasing identification
phaseblock <- unique(BCfile$phaseblock)

# Set the column names
colnames(BCfile) <- c("chr_pos_ref_alt", "ref_count", "alt_count", "tot_depth", "zygosity", "phaseblock", "barcodes")

# These are all of the Barcodes from both the reference and the alternative alleles
barcodes_all <- as.list(strsplit(gsub(";|,", " ", as.character(BCfile$barcodes)), '\\s+'))

Ref_allele_BX <- list()
Alt_allele_BX <- list()

for (index in 1:length(BCfile$chr_pos_ref_alt)) {
  Ref_BX <- unlist(as.list(strsplit(x = BCfile$barcodes[index], split = ",")[[1]])[1])
  if(is.null(Ref_BX)){Ref_BX <- "."}
  
  Alt_BX <- unlist(as.list(strsplit(x = BCfile$barcodes[index], split = ",")[[1]])[2]) 
  if(is.null(Alt_BX)){Alt_BX <- "."}
  
  Ref_allele_BX <- rbind(Ref_allele_BX, Ref_BX)
  Alt_allele_BX <- rbind(Alt_allele_BX, Alt_BX)
}

# Barcode list for the reference allele and the alternative allele
BCfile$Ref_Barcodes <- Ref_allele_BX
BCfile$Alt_Barcodes <- Alt_allele_BX

# Create an empty dataframe
outputDF <- data.frame()

# Arching goal for each variant get the count of matching barcodes
for (zndex in 1:length(BCfile$chr_pos_ref_alt)) {
  # Variant of interest
  row_variant = BCfile$chr_pos_ref_alt[zndex]
  variantOfInterest <- subset(BCfile, grepl((row_variant), BCfile$chr_pos_ref_alt))
  
  # Count of support
  voi_ref_count <- as.integer(variantOfInterest$ref_count)
  voi_alt_count <- as.integer(variantOfInterest$alt_count)
  
  # Alternative allele bar codes for variant of interest
  VOI_altAllele_BC <- unlist(variantOfInterest$Alt_Barcodes)
  VOI_altAllele_BC <- unlist(strsplit(gsub(";", " ", as.character(VOI_altAllele_BC)), '\\s+'))
  VOI_altAllele_BC <- gsub("\\-.*", "", VOI_altAllele_BC)
  Number_of_Barcodes_VOI_alt <- length(VOI_altAllele_BC)
  
  # Reference allele bar codes for variant of interest
  VOI_RefAllele_BC <- unlist(variantOfInterest$Ref_Barcodes)
  VOI_RefAllele_BC <- unlist(strsplit(gsub(";", " ", as.character(VOI_RefAllele_BC)), '\\s+'))
  VOI_RefAllele_BC <- gsub("\\-.*", "", VOI_RefAllele_BC)
  
  Number_of_Barcodes_VOI_Ref <- length(VOI_RefAllele_BC)
  
  # Variant for pairwise comparisons
  comparisonVariants <- subset(BCfile, !grepl((row_variant), BCfile$chr_pos_ref_alt))
  
  # Go through the comparison between the VOI and each of CV
  for (xindex in 1:length(comparisonVariants$chr_pos_ref_alt)) {
    CV_identity <- comparisonVariants$chr_pos_ref_alt[xindex]
    CV_alt_count <- comparisonVariants$alt_count[xindex]
    CV_ref_count <- as.integer(comparisonVariants$ref_count[xindex])
    
    CV_altAllele_BC <- unlist(comparisonVariants$Alt_Barcodes[xindex])
    CV_altAllele_BC <- gsub("\\-.*", "", unlist(strsplit(gsub(";", " ", as.character(CV_altAllele_BC)), '\\s+')))
    
    CV_refAllele_BC <- unlist(comparisonVariants$Ref_Barcodes[xindex])
    CV_refAllele_BC <- gsub("\\-.*", "", unlist(strsplit(gsub(";", " ", as.character(CV_refAllele_BC)), '\\s+')))
    
    # Total barcodes for the variant of interest
    Total_VOI_BC <- length(c(VOI_altAllele_BC, VOI_RefAllele_BC))
    
    # All possible barcodes between two variants 
    All_possible_BC <- length(unique(c(VOI_altAllele_BC, VOI_RefAllele_BC, CV_altAllele_BC, CV_refAllele_BC)))
    
    # This is the In phase (concordance/discordance)
    INPHASE_BC_alt <- length(intersect(VOI_altAllele_BC, CV_altAllele_BC))
    INPHASE_BC_ref <- length(intersect(VOI_RefAllele_BC, CV_refAllele_BC))
    total_inphase <- INPHASE_BC_alt + INPHASE_BC_ref
    
    # This is the out of phase count (trans)
    OutPhase_BC_VOI_alt <- length(intersect(VOI_altAllele_BC, CV_refAllele_BC))
    OutPhase_BC_VOI_ref <- length(intersect(VOI_RefAllele_BC, CV_altAllele_BC))
    total_outphase <- OutPhase_BC_VOI_alt + OutPhase_BC_VOI_ref
    
    indexVariant_informatio <- list(c(row_variant, voi_ref_count, voi_alt_count, CV_identity, CV_ref_count, CV_alt_count, INPHASE_BC_alt, INPHASE_BC_ref, OutPhase_BC_VOI_alt, OutPhase_BC_VOI_ref, All_possible_BC))
    
    outputDF <- rbind(outputDF, do.call(rbind, indexVariant_informatio))
  }
}

colnames(outputDF) <- c("VOI", "voi_ref_count", "voi_alt_count", "CV", "CV_ref_count", "CV_alt_count", "inphase_alt", "inphase_ref", "outphase_VOI_alt", "outphase_VOI_ref", "All_possible_BC")

outputDF <- na.omit(outputDF)

# Write the file out
write.table(outputDF, file = paste(input[1], "phasing_information.txt", sep = "_"), sep = "\t", quote = FALSE)
