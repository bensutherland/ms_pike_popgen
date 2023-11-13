# Analyze the number of fixed differences between various populations in the northern pike dataset
# B. Sutherland (2023-11-09)

# Clear space
#rm(list=ls())

# Load packages
#install.packages("vcfR")
library("vcfR")

# Set working directory
setwd("~/Documents/00_sbio/UVic_northern_pike/02_fixed_diffs")

# Read in data
vcf <- read.vcfR(file = "../01_polymorph_enrich/Eluc.variants.GATK.iteration.2.b.Full.SNP.GATK_HF_removed.minQ20.mmdp10.mxmdp60.mmc10.mac1.homsumfilt.new.samp.names.recode.vcf")

# Extract genotypes
gt.df <- extract.gt(x = vcf, element = "GT", as.numeric = TRUE) # note: also possible to 'return.alleles' or "0/1" format
gt.df[1:5,1:5]

# View all the sample names
colnames(gt.df)
dim(gt.df)

# Make separate dataframes for each of the target pops
gt_HCK.df <- gt.df[, grep(x = colnames(gt.df), pattern = "NJ")]
gt_SLA.df <- gt.df[, grep(x = colnames(gt.df), pattern = "-S")]
gt_CHT.df <- gt.df[, grep(x = colnames(gt.df), pattern = "CR")]
gt_HOO.df <- gt.df[, grep(x = colnames(gt.df), pattern = "-YR")]

#gt_sel.df <- gt.df[, grep(x = colnames(gt.df), pattern = "NJ|S|CR|YR")]

# What are the samples in each pop?
colnames(gt_HCK.df)
colnames(gt_SLA.df)
colnames(gt_CHT.df)
colnames(gt_HOO.df)

# How many in each pop?
ncol(gt_HCK.df) #  9
ncol(gt_SLA.df) # 11
ncol(gt_CHT.df) # 10
ncol(gt_HOO.df) #  5

# Set up list
genos.list <- list()
genos.list[["HCK"]] <- gt_HCK.df
genos.list[["SLA"]] <- gt_SLA.df
genos.list[["CHT"]] <- gt_CHT.df
genos.list[["HOO"]] <- gt_HOO.df

# Run a loop
# Set the number of records
#num.records <- 100000

genos.df <- NULL; pop.target <- NULL; result.list <- list()
for(i in 1:length(genos.list)){
  
  # Select the target population
  pop.target <- names(genos.list)[i]
  
  print(paste0("Working on population: ", pop.target))
  
  # Extract the geno for target pop
  genos.df <- genos.list[[pop.target]]
  
  # Make a matrix to collect information
  result.df <- matrix(data = NA, nrow = nrow(genos.df), ncol = 5)
  colnames(result.df) <- c("locus.name", "num.homo.ref", "num.het", "num.homo.alt", "num.missing")
  result.df <- as.data.frame(result.df)
  head(result.df)
  
  # Debugging
  #result.df <- result.df[1:num.records, ]
  
  # Calculate number of each per line
  for(l in 1:nrow(result.df)){
    
    result.df[l, "locus.name"]   <- rownames(genos.df)[l]
    result.df[l, "num.homo.ref"] <- sum(genos.df[l, ]=="0", na.rm = T)
    result.df[l, "num.het"]      <- sum(genos.df[l, ]=="1", na.rm = T)
    result.df[l, "num.homo.alt"] <- sum(genos.df[l, ]=="2", na.rm = T)
    result.df[l, "num.missing"]  <- sum(is.na(genos.df[l,]))
    
  }
  
  head(result.df, n = 20)
  
  dim(result.df)
  
  # Check for NA
  if(sum(is.na(c(result.df$num.homo.ref, result.df$num.het, result.df$num.homo.alt)))==0){
    
    print("All looks good, continue")
    
  }else{
    
    print("Something went wrong and a sum is NA, check your formula")
    
  }
  
  # Inspect results
  dim(result.df)
  head(result.df, n = 20)
  
  # Any homozygous alternate? 
  result.df[which(result.df$num.homo.alt!=0), ] # yes, but not many
  
  # Do you have any hets? 
  head(result.df[which(result.df$num.het!=0), ], n = 20)  # yes, many
  
  # Remove heterozygous, as they are not fixed
  result.df <- result.df[result.df$num.het==0, ]
  dim(result.df)
  
  # Discard any loci that have both homozygous ref and alt
  result.df <- result.df[result.df$num.homo.ref!=0 & result.df$num.homo.alt==0 | result.df$num.homo.ref==0 & result.df$num.homo.alt!=0, ]
  dim(result.df)
  
  table(result.df$num.homo.ref)
  sum(table(result.df$num.homo.ref)) # May not always work
  table(result.df$num.homo.alt)
  table(result.df$num.homo.alt)
  
  # result.list <- list()
  # result.list[["HCK"]] <- result.df
  # head(result.list[["HCK"]])
  
  result.list[[pop.target]] <- result.df
  
}

# All the results are now in result.list
names(result.list)

CHT.df <- result.list[["CHT"]]
HOO.df <- result.list[["HOO"]]
HCK.df <- result.list[["HCK"]]
SLA.df <- result.list[["SLA"]]


# Reporting
print("Now only have fixed alleles per pop")
print("Remember, this is considering fixed differences only")


# Inspect results manually
head(CHT.df)
dim(CHT.df)                # 85,267 records
table(CHT.df$num.homo.ref)
table(CHT.df$num.homo.alt) # 2 alt homozyg

head(HOO.df)
dim(HOO.df)                # 86,640 records
table(HOO.df$num.homo.ref)
table(HOO.df$num.homo.alt) # 2 alt homozyg

head(HCK.df)
dim(HCK.df)                # 96,864 records
table(HCK.df$num.homo.ref)
table(HCK.df$num.homo.alt) # 0 alt homozyg

head(SLA.df)
dim(SLA.df)                # 98,154 records
table(SLA.df$num.homo.ref)
table(SLA.df$num.homo.alt) # 0 alt homozyg


#### Contrasts, fixed diffs ####
# Set the variables for pop 1 or pop 2, then produce a report below
pop1 <- "CHT"
pop2 <- "SLA"

# How many variants are fixed for the ref or alt allele in pop 1?
sum(result.list[[pop1]]$num.homo.ref > 0) # 85,265
sum(result.list[[pop1]]$num.homo.alt > 0) #      2

# How many variants are fixed for the ref allele in pop 2?
sum(result.list[[pop2]]$num.homo.ref > 0) # 98,154
sum(result.list[[pop2]]$num.homo.alt > 0) #      0

# and what are they? 
fixed_ref_loci_pop1.vec <- result.list[[pop1]][result.list[[pop1]]$num.homo.ref > 0 , "locus.name"]
fixed_alt_loci_pop1.vec <- result.list[[pop1]][result.list[[pop1]]$num.homo.alt > 0 , "locus.name"]

# and what are they? 
fixed_ref_loci_pop2.vec <- result.list[[pop1]][result.list[[pop2]]$num.homo.ref > 0 , "locus.name"]
fixed_alt_loci_pop2.vec <- result.list[[pop1]][result.list[[pop2]]$num.homo.alt > 0 , "locus.name"]

# How many of the fixed alt in pop 1 are fixed ref in pop 2? 
length(intersect(x = fixed_alt_loci_pop1.vec, y = fixed_ref_loci_pop2.vec)) # 2

# How many of the fixed alt in pop 2 are fixed ref in pop 1? 
length(intersect(x = fixed_alt_loci_pop2.vec, y = fixed_ref_loci_pop1.vec)) # 0

# With these items, it should provide all of the information you need for the study
