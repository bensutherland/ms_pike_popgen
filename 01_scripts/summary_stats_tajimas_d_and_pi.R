# Calculate summary stats for Tajima's D and nucleotide diversity (pi) per population
# B. Sutherland (2023-11-27)

# Clear space
#rm(list=ls())

# Load packages

# Set working directory to the ms_scallop_popgen repo
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

#### 01. Tajima's D ####
# List input files
input_files.vec <- list.files(path = "03_results/", pattern = "Tajima.D")

df <- NULL
for(i in 1:length(input_files.vec)){
  
  # Reporting
  print(paste0("Average Tajima's D for file: ", input_files.vec[i]))
    
  # Read in data
  df  <-  read.delim(file = paste0("03_results/", input_files.vec[i]))
  
  # Calculate mean
  print("**Average: ")
  print(mean(df$TajimaD, na.rm = T))
  
}


#### 02. Nucleotide Diversity (pi) ####
# List input files
input_files.vec <- list.files(path = "03_results/", pattern = ".sites.pi")

df <- NULL
for(i in 1:length(input_files.vec)){
  
  # Reporting
  print(paste0("Average sites-pi for file: ", input_files.vec[i]))
  
  # Read in data
  df  <-  read.delim(file = paste0("03_results/", input_files.vec[i]))
  
  # Calculate mean
  print("**Average: ")
  print(mean(df$PI, na.rm = T))
  
}
