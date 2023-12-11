# Per population stats
# Ben J. G. Sutherland, 2023-11-23

# Clear space
#rm(list=ls())

# Load packages
#install.packages("dartR")
#install.packages("vcfR")
#install.packages("BiocManager")
#install.packages("poppr")
library("BiocManager")
#BiocManager::install("SNPRelate")
library("SNPRelate")
library("dartR")
library("vcfR")
library("poppr")

# Set working directory to the ms_scallop_popgen repo
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# Set variables
options(scipen = 9999999)
input_VCF.FN <- "02_input_data/Eluc.variants.GATK.iteration.2.b.Full.SNP.GATK_HF_removed.minQ20.mmdp10.mxmdp60.mmc10.mac1.homsumfilt.new.samp.names.recode_biallelic_only.vcf"


#### 00. Read in data ####
# Read in VCF
my_data.vcf <- read.vcfR(file = input_VCF.FN)
my_data.vcf

#head(my_data.vcf@fix)
#my_data.vcf@fix


#### 01. Prepare different formats ####
# Convert VCF to genotype dataframe
gt.df  <- extract.gt(x = my_data.vcf, element = "GT")
gt.df[1:5,1:5]
head(rownames(gt.df))

# # Optional - check to confirm that multiallelic variants are gone: 
# for(i in 1:ncol(gt.df)){
#   
#   print(i)
#   print(table(gt.df[,i], useNA = "ifany"))
#   
# }
## Looks good

# Convert VCF to genind
my_data.gid <- vcfR2genind(x = my_data.vcf)
my_data.gid
head(locNames(my_data.gid))

# Compare locus names bw geno df and genind
length(intersect(x =  locNames(my_data.gid), y = (rownames(gt.df))))
nLoc(my_data.gid)

# Assign population attribute of genind based on individual names
inds <- indNames(my_data.gid)
inds[grep(pattern = "-CL-", x = inds)] <- "CHA"
inds[grep(pattern = "-CaG-", x = inds)] <- "CAS"
inds[grep(pattern = "-Mb", x = inds)] <- "WHI"
inds[grep(pattern = "-NJ", x = inds)] <- "HCK"
inds[grep(pattern = "-S", x = inds)] <- "SLA"
inds[grep(pattern = "-PL", x = inds)] <- "PAL"
inds[grep(pattern = "-CR", x = inds)] <- "CHT"
inds[grep(pattern = "-YR", x = inds)] <- "HOO"
inds <- as.character(inds)
inds

pop(my_data.gid) <- inds
rm(inds)

# Confirm OK: 
as.data.frame(pop(my_data.gid), indNames(my_data.gid))
table(pop(my_data.gid))

# Convert genind to genlight
my_data.gl <- gi2gl(gi = my_data.gid, parallel = T, verbose = T)
my_data.gl


#### 02. Per population mean HOBS and FIS ####
# HOBS and FIS
df <- gl.report.heterozygosity(x = my_data.gl, method = "pop")
write.table(x = df, file = "03_results/HOBS_summary_stats.txt", quote = F
            , sep = "\t", row.names = F)


#### 03. Private alleles ####
# Calculate private alleles per population
pa <- private_alleles(gid = my_data.gid) # calc'd by pop by default
dim(pa) 
pa[1:4,1:10]    # rows: population; cols: private alleles

# Transpose pa output
pa.t <- t(pa)   # Summary table of pa per pop (obs number of pa to each pop)
dim(pa.t)
pa.t[1:10,]     # rows: private alleles; cols = pop
# note: the output is a unit per allele, e.g., a single het in the pop will be val of 1, two hets = 2, two homo alt = 4

# To confirm the above, inspect the first row of private allele output in the genotype df (note: drop allelic designation)
gt.df["LG07_3041", ] # there are two hets in CHT
gt.df["LG07_5166", ] # there are two hets in HCK

# Remove allele ID from the pa locus name
rownames(x = pa.t) <- gsub(pattern = "\\.[0-9]", x = rownames(x = pa.t), replacement = "", perl = T)
head(pa.t)

# Data check: any duplicates after removing allele? (there shouldn't be)
table(duplicated(rownames(pa.t))) # all false, OK

# Identify how many private alleles per pop
for(i in 1:ncol(pa.t)){
  
  print(paste0("pop: ", colnames(pa.t)[i]))
  
  print(table(pa.t[,i] > 0))
  
}

# write.csv(x = pa, file = "03_results/private_alleles.csv", quote = F)

# If want to limit to pops with at least 5 indiv
table(pop(my_data.gid))

# What are the tallies of pa per pop for those with at least five indiv? 
head(pa.t)
table(pa.t[,"SLA"])
table(pa.t[,"CHT"])
table(pa.t[,"HCK"])
table(pa.t[,"WHI"])
table(pa.t[,"HOO"])

# Make df of all pa markers, and the name of their pop
head(pa.t)
result.df <- matrix(data = "NA", nrow = nrow(pa.t), ncol = 2)
result.df <- as.data.frame(result.df)
colnames(result.df) <- c("locus_name", "pa_pop")

# What variants are pa in each pop? 
head(pa.t)

pa.list <- list()
for(i in 1:ncol(pa.t)){
  
  pa.list[[colnames(pa.t)[i]]] <- rownames(pa.t)[pa.t[,i] > 0]
  
}

pa.df <- NULL; slice <- NULL; poi <- NULL
for(i in 1:length(pa.list)){
  
  poi <- names(pa.list)[i]
  slice <- pa.list[[i]]
  
  slice <- cbind(poi, slice)
  slice <- as.data.frame(slice)
  colnames(slice) <- c("pop", "pa.locus")
  
  pa.df <- rbind(slice, pa.df)
  
}

head(pa.df)
table(pa.df$pop)
# This now contains a pop column, and the locus name for the private allele
# this can be compared to the fixed allele information from 'fixed_diffs*.R' 

write.table(x = pa.df, file = "03_results/private_alleles_with_pop.txt", quote = F, sep = "\t"
            , row.names = F)

### NEXT: ###
# - how many of the private alleles are fixed? This may use 'fixed diffs' script
# - compare directly to the number of variable loci per pop

