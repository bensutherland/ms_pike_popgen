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
# Read in VCF as a genlight file
input.gl <- gl.read.vcf(vcffile = input_VCF.FN)

my_data.vcf <- read.vcfR(file = input_VCF.FN)
my_data.vcf

# Extract genotypes
gt.df <- extract.gt(x = my_data.vcf, element = "GT"
) 
gt.df[1:5,1:5]

# multiallelic gone? 
for(i in 1:ncol(gt.df)){
  
  print(i)
  print(table(gt.df[,i], useNA = "ifany"))
  
}

# Looks good

# Convert the VCF to genind
my_data.gid <- vcfR2genind(x = my_data.vcf)

# Assign population
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


pop(my_data.gid) <- inds
rm(inds)

# Confirm OK: 
as.data.frame(pop(my_data.gid), indNames(my_data.gid))


table(pop(my_data.gid))

# Convert genind to genlight
my_data.gl <- gi2gl(gi = my_data.gid, parallel = T, verbose = T)
my_data.gl


# HOBS and FIS
df <- gl.report.heterozygosity(x = my_data.gl, method = "pop")
write.table(x = df, file = "03_results/HOBS_summary_stats.txt", quote = F
            , sep = "\t", row.names = F)



#### Private alleles  ####
pa <- private_alleles(gid = my_data.gid)
pa.t <- t(pa) # Summary table of pa per pop (obs number of pa to each pop)
# note that the number will give one per allele, e.g., a single het in the pop will be val of 1, two hets = 2, two homo alt = 4
gt.df["LG07_3041", ] # to demonstrate

head(pa.t)
for(i in 1:ncol(pa.t)){
  
  print(paste0("pop: ", colnames(pa.t)[i]))
  
  print(table(pa.t[,i] > 0))
  
}
table(pa.t[,"CHA"] > 0)  #   T: 14,736
write.csv(x = pa, file = "03_results/private_alleles.csv", quote = F)

### NEXT: ###
# - limit to only populations with a specific number of individuals
# - how many of the private alleles are fixed? This may use 'fixed diffs' script
# - compare directly to the number of variable loci per pop (this may be avail already)
