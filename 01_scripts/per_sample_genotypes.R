# Summarize proportions or counts of each genotype per individual
# B. Sutherland (2023-11-14)

# Clear space
#rm(list=ls())

# Load packages
#install.packages("vcfR")
#install.packages("rstudioapi")
#install.packages("dplyr")
#install.packages("tidyr")
library("tidyr")
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
head(summary_all.df)

# Restrict to only the target genos for plotting
summary_all.df <- summary_all.df[which(summary_all.df$geno=="0/1" | 
                                         summary_all.df$geno=="1/1"), ]

# Rename samples
summary_all.df$indiv

summary_all.df$indiv <- gsub(pattern = "Eluc", replacement = "", x = summary_all.df$indiv) # drop Eluc
summary_all.df$indiv <- gsub(pattern = "Elub", replacement = "", x = summary_all.df$indiv) # typo in data
summary_all.df$indiv
summary_all.df$indiv <- gsub(pattern = "-YR", replacement = "HOO_", x = summary_all.df$indiv)
summary_all.df$indiv <- gsub(pattern = "-Mb", replacement = "WHI_", x = summary_all.df$indiv)
summary_all.df$indiv <- gsub(pattern = "-CR", replacement = "CHT_", x = summary_all.df$indiv)
summary_all.df$indiv <- gsub(pattern =  "-S", replacement = "SLA_", x = summary_all.df$indiv)
summary_all.df$indiv <- gsub(pattern = "-PL", replacement = "PAL_", x = summary_all.df$indiv)
summary_all.df$indiv <- gsub(pattern = "-NJ", replacement = "HCK_", x = summary_all.df$indiv)
summary_all.df$indiv <- gsub(pattern = "-CL", replacement = "CHA_1", x = summary_all.df$indiv) # note: added 1 for matching
summary_all.df$indiv <- gsub(pattern = "-CaG", replacement = "CAS_1", x = summary_all.df$indiv) # note: added 1 for matching
#summary_all.df$indiv <- gsub(pattern = "-", replacement = "_", x = summary_all.df$indiv)
summary_all.df$indiv <- gsub(pattern = "-.*", replacement = "", x = summary_all.df$indiv)

summary_all.df$indiv

# Reorder
summary_all.df$order <- NA
summary_all.df$order[grep(pattern = "CHT", x = summary_all.df$indiv)] <- 1
summary_all.df$order[grep(pattern = "HOO", x = summary_all.df$indiv)] <- 2
summary_all.df$order[grep(pattern = "PAL", x = summary_all.df$indiv)] <- 3
summary_all.df$order[grep(pattern = "CHA", x = summary_all.df$indiv)] <- 4
summary_all.df$order[grep(pattern = "CAS", x = summary_all.df$indiv)] <- 5
summary_all.df$order[grep(pattern = "WHI", x = summary_all.df$indiv)] <- 6
summary_all.df$order[grep(pattern = "SLA", x = summary_all.df$indiv)] <- 7
summary_all.df$order[grep(pattern = "HCK", x = summary_all.df$indiv)] <- 8

# Can't figure out how to use backreference to 
# grep(pattern = '_[0-9]$', x = summary_all.df$indiv, perl = T, invert = T)
# #gsub(pattern = '_', replacement = "_0", x = summary_all.df$indiv[grep(pattern = '_[0-9]$', x = summary_all.df$indiv, perl = T)])
# gsub(pattern = '_[0-9]$', x = summary_all.df$indiv, replacement = "_0", fixed = F)

# Need to reorder by pop
# TODO: update axis labels

#summary_all.df.bck <- summary_all.df
# summary_all.df <- summary_all.df[which(summary_all.df$geno=="0/0" | summary_all.df$geno=="0/1" | 
#                                 summary_all.df$geno=="1/1"), ]

# Add leading zero
summary_all.df <- separate(data = summary_all.df, col = "indiv", into = c("pop", "num"), sep = "_", remove = T)
summary_all.df$num <- as.numeric(summary_all.df$num)
summary_all.df$num <- sprintf(fmt = '%02d', summary_all.df$num)
summary_all.df$indiv <- paste0(summary_all.df$pop, "_", summary_all.df$num)
head(summary_all.df)
summary_all.df <- summary_all.df[,c("indiv", "geno", "count", "order")]

# # Rename genos (would implement here)
# summary_all.df$geno <- gsub(pattern = "0/1", replacement = "Het.", x = summary_all.df$geno)

# Plot
options(scipen = 99999999)

#install.packages("ggplot2")
library("ggplot2")

# Reorder attempt
p <- ggplot(data = summary_all.df, aes(fill=geno, x = reorder(indiv, order), y = count)) +
       geom_bar(position='stack', stat = 'identity') +
       theme_bw() +
       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
             #, panel.border = element_blank()
             , panel.grid.major = element_blank()
             , panel.grid.minor = element_blank()
             #,
             ) +
       scale_fill_grey(start = 0.3, end = 0.7) + # grey colour
       xlab("Individual") + 
       ylab("Number variants") 
p

# Save plot
pdf(file = "03_results/num_var_per_indiv.pdf", width = 9, height = 3.5)
print(p)
dev.off()

variants_per_indiv.plot <- p

# Save per sample genotype count plot as object
save(variants_per_indiv.plot, file = "03_results/per_ind_genos.Rdata")

# Calculate per population summary stats
head(summary_all.df)

sumstats.df  <- separate(data = summary_all.df, col = "indiv", into = c("pop", "indiv")
                       , sep = "_", remove = F)

head(sumstats.df)
sumstats.df$indname <- paste0(sumstats.df$pop, "_", sumstats.df$indiv)
head(sumstats.df)

# Separate the genos into two df then recombine
het.df  <- sumstats.df[sumstats.df$geno=="0/1", ]
homo.df <- sumstats.df[sumstats.df$geno=="1/1", ]
head(het.df)
head(homo.df)

wide.df <- merge(x = het.df, y = homo.df, by = "indname", all = T)
head(wide.df)
wide.df <- wide.df[,c("pop.x", "indname", "geno.x", "count.x", "geno.y", "count.y")]
head(wide.df)
colnames(wide.df) <- c("pop", "ind", "geno.het", "count.het", "geno.homo.alt", "count.homo.alt")
head(wide.df)
write.table(x = wide.df, file = "03_results/number_genos_per_indiv.txt", quote = F, sep = "\t", row.names = F)


# Calculate the average number of each genotype per individual within each population
mean(sumstats.df[sumstats.df$geno=="0/1" & sumstats.df$pop=="HOO", "count"])
aggregate(count ~ pop, data = sumstats.df[sumstats.df$geno=="0/1",], FUN = mean)
aggregate(count ~ pop, data = sumstats.df[sumstats.df$geno=="1/1",], FUN = mean)

# Complete
