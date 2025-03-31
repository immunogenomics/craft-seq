#!/usr/bin/env Rscript
# Date: 02/06/2025
# Author: Michelle Curtis, Yuriy Baglaenko
# This script converts individual allele tables resulting from running CRISPREsso2 on each cell into a single allele table matrix.

source("../misc_funcs.R")
library(stringr)

DNAfiles = "./DNAOutput_R2/"

#Alleles Function for Import - Updated- Fixed Naming to Automate and Changed to Capture the top 10 alleles per cell. 
AllelesFunction <- 
function (filepath, number_of_alleles = 10) #Filepath is the folder with all the allele tables. 
{
    filepath <- sub("/$", "", filepath)
    files <- list.files(path = filepath, pattern = "Alleles_frequency_table.txt$", full.names = T, 
        recursive = T) 
        
    matrix <- rbindlist(lapply(files, function(f) {
    # Extract Plate_Well ID from the file path
    plate_id <- str_extract(f, "[^/]+(?=/CRISPRessoBatch_on_batch/CRISPResso_on_[^/]+)") %>%
                paste0("_", str_extract(f, "(?<=CRISPResso_on_)[^/]+"))
    
        # Read file and add Plate_Well column
        mat_file <- fread(f, nrows = number_of_alleles, select = c(1, 5, 6, 7, 8, 9))
        mat_file[, Plate_Well := plate_id]  # Assign the extracted Plate_Well ID

        return(mat_file)  # Return the modified table
    }), use.names = TRUE, fill = TRUE)  # Bind all tables together

    matrix = cbind(matrix[, (ncol(matrix)):(ncol(matrix))], matrix[, 1:(ncol(matrix)-1)])
    
    return(matrix)
}

#Upload all Alleles
alleles_DNA <- AllelesFunction(DNAfiles)
alleles_DNA <- alleles_DNA %>% separate(Plate_Well, into = c("DNA_Barcode", "Well_ID"))

write_csv(alleles_DNA, "alleles_DNA.csv")