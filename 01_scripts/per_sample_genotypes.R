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
gt.df <- extract.gt(x = vcf, element = "GT", as.numeric = TRUE
                    ) # note: also possible to 'return.alleles' or "0/1" format
gt.df[1:5,1:5]

# View all the sample names
colnames(gt.df)
dim(gt.df)

# Manual inspection
table(gt.df[,"Eluc-PL1-U"])
#table(gt.df[,"Eluc-CaG-F"])

test.df <- extract.gt(x = vcf, element = "GT")
head(test.df)
table(test.df[,"Eluc-PL1-U"], useNA = "ifany")


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
  #print(i)
  
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

# Plot
options(scipen = 99999999)

#install.packages("ggplot2")
library("ggplot2")

ggplot(data = summary_all.df, aes(fill=Var1, y = Freq, x = name.vec)) + 
  geom_bar(position='stack', stat = 'identity')


#### 02. Analyze ####
# Characterize features per pop in loop
genos.df <- NULL; pop.target <- NULL; result.list <- list()
for(i in 1:ncol(gt.df)){
  
  # Select the target population
  pop.target <- names(genos.list)[i]
  
  print(paste0("Working on population: ", pop.target))
  
  # Extract the geno for the target pop
  genos.df <- genos.list[[pop.target]]
  
  # Prepare blank matrix to collect information
  result.df <- matrix(data = NA, nrow = nrow(genos.df), ncol = 5)
  colnames(result.df) <- c("locus.name", "num.homo.ref", "num.het", "num.homo.alt", "num.missing")
  result.df <- as.data.frame(result.df)
  head(result.df)
  
  # Characterize by row
  ### OLD METHOD ###
  # for(l in 1:nrow(result.df)){
  #   
  #   result.df[l, "locus.name"]   <- rownames(genos.df)[l]
  #   result.df[l, "num.homo.ref"] <- sum(genos.df[l, ]=="0", na.rm = T)
  #   result.df[l, "num.het"]      <- sum(genos.df[l, ]=="1", na.rm = T)
  #   result.df[l, "num.homo.alt"] <- sum(genos.df[l, ]=="2", na.rm = T)
  #   result.df[l, "num.missing"]  <- sum(is.na(genos.df[l,]))
  #   
  # }
  ### /END/ OLD METHOD ###
  
  result.df$locus.name   <- rownames(genos.df) # locus name
  result.df$num.homo.ref <- rowSums(genos.df=="0", na.rm = T) # homo ref
  result.df$num.het      <- rowSums(genos.df=="1", na.rm = T) # het
  result.df$num.homo.alt <- rowSums(genos.df=="2", na.rm = T) # homo alt
  result.df$num.missing  <- rowSums(is.na(genos.df))          # NA
  
  #head(result.df, n = 20)
  #dim(result.df)
  
  # Data checking
  if(sum(is.na(c(result.df$num.homo.ref, result.df$num.het, result.df$num.homo.alt)))==0){
    
    print("All looks good, continue")
    
  }else{
    
    print("Something went wrong and a sum is NA, check your formula")
    
  }
  
  print(paste0("This population has ", nrow(result.df), " variants before any are removed."))
  
  # Any homozygous alternate? 
  #result.df[which(result.df$num.homo.alt!=0), ]
  
  # Do you have any hets? 
  #result.df[which(result.df$num.het!=0), ]
  
  # Only retain variants with exactly zero heterozygous genotypes, as looking for fixed diffs
  print("Removing any variants that have heterozygous genotypes in this population.")
  result.df <- result.df[result.df$num.het==0, ]
  print(paste0("After removing hets, there are ", nrow(result.df), " variants remaining."))
  
  # Discard any loci that have both homozygous ref and alt
  print("Keeping variants that have homozygous ref genotypes and no homozygous alts, or visa-versa, within this target population")
  result.df <- result.df[result.df$num.homo.ref!=0 & result.df$num.homo.alt==0 | result.df$num.homo.ref==0 & result.df$num.homo.alt!=0, ]
  print(paste0("After removing unfixed variants, there are ", nrow(result.df), " variants remaining."))
  
  # table(result.df$num.homo.ref)
  # sum(table(result.df$num.homo.ref)) # May not always work
  # table(result.df$num.homo.alt)
  # table(result.df$num.homo.alt)
  
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
head(HCK.df)
dim(HCK.df)                # 1,074,768 records
table(HCK.df$num.homo.ref)
table(HCK.df$num.homo.alt) #         4 fixed alt homozyg

head(SLA.df)
dim(SLA.df)                # 1,089,447 records
table(SLA.df$num.homo.ref)
table(SLA.df$num.homo.alt) #         6 fixed alt homozyg

head(CHT.df)
dim(CHT.df)                #   887,207 records
table(CHT.df$num.homo.ref) 
table(CHT.df$num.homo.alt) #        12 alt homozyg

head(HOO.df)
dim(HOO.df)                #   958,025 records
table(HOO.df$num.homo.ref)
table(HOO.df$num.homo.alt) #        14 alt homozyg


#### Contrasts, fixed diffs ####
# Set the variables for pop 1 or pop 2, then produce a report below
# ## main contrast
# pop1 <- "CHT"
# pop2 <- "SLA"

# ## within west
# pop1 <- "CHT"
# pop2 <- "HOO"

## within east
pop1 <- "SLA"
pop2 <- "HCK"


# How many variants are fixed for the ref or alt allele in pop 1?
sum(result.list[[pop1]]$num.homo.ref > 0)
sum(result.list[[pop1]]$num.homo.alt > 0)

# How many variants are fixed for the ref allele in pop 2?
sum(result.list[[pop2]]$num.homo.ref > 0)
sum(result.list[[pop2]]$num.homo.alt > 0)

# and what are they (pop 1)? 
fixed_ref_loci_pop1.vec <- result.list[[pop1]][result.list[[pop1]]$num.homo.ref > 0 , "locus.name"]
fixed_alt_loci_pop1.vec <- result.list[[pop1]][result.list[[pop1]]$num.homo.alt > 0 , "locus.name"]

# and what are they (pop 2)? 
fixed_ref_loci_pop2.vec <- result.list[[pop1]][result.list[[pop2]]$num.homo.ref > 0 , "locus.name"]
fixed_alt_loci_pop2.vec <- result.list[[pop1]][result.list[[pop2]]$num.homo.alt > 0 , "locus.name"]

# How many of the fixed alt in pop 1 are fixed ref in pop 2? 
length(intersect(x = fixed_alt_loci_pop1.vec, y = fixed_ref_loci_pop2.vec))

# How many of the fixed alt in pop 2 are fixed ref in pop 1? 
length(intersect(x = fixed_alt_loci_pop2.vec, y = fixed_ref_loci_pop1.vec))

# With these items, it should provide all of the information you need for the study


