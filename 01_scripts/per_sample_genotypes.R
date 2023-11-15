# Summarize proportions or counts of each genotype per individual
# B. Sutherland (2023-11-14)

# Clear space
#rm(list=ls())

# Load packages
#install.packages("vcfR")
#install.packages("rstudioapi")
#install.packages("dplyr")
library("vcfR")
library("rstudioapi")
library("dplyr")

# Set working directory to the ms_scallop_popgen repo
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

#### 01. Set up ####
# Set filenames
vcf.FN <- "02_input_data/Eluc.variants.GATK.iteration.2.b.Full.SNP.GATK_HF_removed.minQ20.mmdp10.mxmdp60.mmc10.mac1.homsumfilt.new.samp.names.recode.vcf"

# Read in data
vcf <- read.vcfR(file = vcf.FN)

# Extract genotypes
gt.df <- extract.gt(x = vcf, element = "GT"
                    ) # note: do not use the 'as.numeric' method because there are more than two alleles in the VCF and this causes issues
gt.df[1:5,1:5]

# View all the sample names
colnames(gt.df)
dim(gt.df)

# Manual inspection
table(gt.df[,"Eluc-PL1-U"], useNA = "ifany")
#table(gt.df[,"Eluc-CaG-F"], useNA = "ifany")

# Inspect
gt.df[1:5,1:5]
gt.df <- as.data.frame(gt.df)
str(gt.df)

## Proof of principle
# test.df <- as.data.frame(table(gt.df$`Eluc-CL-M`, useNA = "ifany"))
# test.df

# Create summary df
result <- NULL; summary.df <- NULL; summary_all.df <- NULL
for(i in 1:ncol(gt.df)){
  
  # Debugging
  print(i)
  
  # Summarize using table
  result <- as.data.frame(table(gt.df[,i], useNA = "ifany"))
  
  # Create identifier
  name.vec <- rep(colnames(gt.df)[i], times = nrow(result))
  
  # Combine identifier and data
  summary_indiv.df <- cbind(name.vec, result)
  
  # Combine with other samples
  summary_all.df   <- rbind(summary_indiv.df, summary_all.df)
  
}

summary_all.df
summary_all.df <- as.data.frame(summary_all.df)
summary_all.df[1:15,]

colnames(summary_all.df)
colnames(summary_all.df) <- c("indiv", "geno", "count")

# Need to reorder by pop
# TODO: use gsub and grep on a vector of sample names to provide order (do further up)
# TODO: replace namings with the new pop namings
# TODO: specify colours
# TODO: try a version with only 0/1 and 1/1 values
# TODO: remove gray grid
# TODO: update axis labels

#summary_all.df.bck <- summary_all.df
# summary_all.df <- summary_all.df[which(summary_all.df$geno=="0/0" | summary_all.df$geno=="0/1" | 
#                                 summary_all.df$geno=="1/1"), ]

summary_all.df <- summary_all.df[which(summary_all.df$geno=="0/1" | 
                                         summary_all.df$geno=="1/1"), ]


# Plot
options(scipen = 99999999)

#install.packages("ggplot2")
library("ggplot2")

p <- ggplot(data = summary_all.df, aes(fill=geno, y = count, x = indiv)) + 
       geom_bar(position='stack', stat = 'identity') +
       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p

