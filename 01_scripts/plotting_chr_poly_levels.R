# Visualizing the differences in chr polymorphism levels
# Ben J. G. Sutherland, 2023-11-01

# Clear space
rm(list=ls())

# Set working directory and options
setwd("~/Documents/00_sbio/UVic_northern_pike/01_polymorph_enrich")
options(scipen = 9999999)

# Load libraries
library(ggplot2)

# Question: are there significantly enriched polymorphic regions in the genome 
#  (e.g., chr LG09, LG11, LG24) as described in Johnson et al. preprint?


# Read in data
coverage.df <- read.delim(file = "coverage.txt", header = F)
head(coverage.df, n = 10)
colnames(coverage.df) <- c("LG", "start", "stop", "poly.sum")
dim(coverage.df)
str(coverage.df)
coverage.df$start <- as.numeric(coverage.df$start)
coverage.df$stop  <- as.numeric(coverage.df$stop)

# # Keep only the chromosomes (named as 'LG')
# unique(coverage.df$LG)
# coverage.df <- coverage.df[grep(pattern = "LG", x = coverage.df$LG), ]
# unique(coverage.df$LG)

# Sort by chr
coverage.df <- coverage.df[with(coverage.df, order(coverage.df$LG, coverage.df$start)), ]
unique(coverage.df$LG)
dim(coverage.df)
head(coverage.df)

# Generate a plotting variable
#coverage.df$plotting.poly.sum <- coverage.df$poly.sum / 1000 


# Basic barplot of poly sum per chr
p<-ggplot(data=coverage.df, aes(x=LG, y=poly.sum)) +
  geom_bar(stat="identity")
p

head(coverage.df)

#### Polymorphism by chr length ####
# Determine the length of each chr
len.chr = NULL; len.chr.list <- list()
for(i in 1:length(unique(coverage.df$LG))){

  chr.o.i <- unique(coverage.df$LG)[i]
  print(chr.o.i)

  len.chr.list[[chr.o.i]] <- max(coverage.df[coverage.df$LG==chr.o.i, "stop"])

}

len.chr.list
len.chr.df <- as.data.frame(unlist(len.chr.list))
len.chr.df$chr <- rownames(len.chr.df)
colnames(len.chr.df) <- c("len", "chr")
len.chr.df <- len.chr.df[,c("chr", "len")]
#len.chr.df$cumsum <- cumsum(len.chr.df$len)
head(len.chr.df)
tail(len.chr.df)

# Attach the polymorphism level to the chromosomes
coverage.df$LG
poly_by_chr.df <- aggregate(coverage.df$poly.sum, by=list(LG=coverage.df$LG), FUN=sum)
str(poly_by_chr.df)
head(poly_by_chr.df)

# Confirm ok? 
sum(coverage.df[coverage.df$LG=="LG01", "poly.sum"]) # OK!

head(len.chr.df)
head(poly_by_chr.df)
chr_length_and_cumul_poly.df <- merge(x = len.chr.df, y = poly_by_chr.df, by.x = "chr", by.y = "LG")
head(chr_length_and_cumul_poly.df)
colnames(chr_length_and_cumul_poly.df) <- c("chr", "length", "poly.sum")
head(chr_length_and_cumul_poly.df)
chr_length_and_cumul_poly.df$length_per.kb <- chr_length_and_cumul_poly.df$length / 1000

plot(chr_length_and_cumul_poly.df$poly.sum / chr_length_and_cumul_poly.df$length_per.kb
     , xaxt = "n"
     , ylab = "avg. SNP / kb"
     , las = 1
     , xlab = "chr"
     )
axis(side = 1, at = seq(1:nrow(chr_length_and_cumul_poly.df)))



#### Approach to calculate number of bins ####
bins <- table(coverage.df$LG)
bins.df <- as.data.frame(names(bins))
colnames(bins.df) <- "chr"
bins.df$bins <- as.numeric(bins)
bins.df$cumsum <- cumsum(bins.df$bins)

head(bins.df)

bins.df$start.bin <- NA; bins.df$end.bin <- NA; bins.df$label.pos <- NA
for(i in 1:nrow(bins.df)){
  
  if(i == 1){
    
    bins.df$start.bin[i] <- 0
    
  }else{
    
    bins.df$start.bin[i] <- bins.df$cumsum[i-1]
    
  }
  
  bins.df$end.bin[i] <- bins.df$cumsum[i]
  
  bins.df$label.pos[i] <- ((bins.df$end.bin[i] - bins.df$start.bin[i])/2) + bins.df$start.bin[i]

}

head(bins.df)

# Create colour vector for plotting
bins.df$plot.colour <- "black"
bins.df$plot.colour[c(FALSE, TRUE)] <- "darkgrey"
bins.df


### TEST ### ### NOT WORKING ###
#head(bins.df)
# xlim <- c(0, length(bins.df$end.bin)*1.25)
#xlim <- c(0, nrow(coverage.df)*1.25)

#barplot(coverage.df$poly.sum, xlim = xlim)

head(coverage.df)
head(bins.df)

# TODO: add colours
coverage.df$plot.color <- NA
for(i in 1:length(coverage.df$plot.color)){
  
  coverage.df$plot.color[i] <- bins.df[bins.df$chr %in% coverage.df$LG[i], "plot.colour"]
  
}

head(coverage.df, n = 10)

pdf(file = "plot_SNP_sum_per_window.pdf", width = 15, height = 4)
plot(x = seq(1:nrow(coverage.df)), y = coverage.df$poly.sum
     , xaxt = "n"
     , xlab = "Chromosome"
     , ylab = "Per window SNP sum"
     , las = 1
     , pch = 16
     , col = coverage.df$plot.color
     , cex = 1
     )
abline(v = bins.df$end.bin, lty = 2)
abline(v = bins.df$start.bin, lty = 2)
axis(side = 1, at = bins.df$label.pos, labels = bins.df$chr
     , tick = T, las = 2
     )

abline(h = min(outlier_polymorphism.vec))

dev.off()

# #### Plot ####
# # Plot an empty plot
# #par(mfrow=c(1,1), mar= c(2,3,0.5,1) + 0.2, mgp = c(2,0.75,0))
# plot(x = c(0, max(cumul.leng)), y = c(0, max(gwas.fst$fst+0.025)), type = "n"
#      , xaxt = 'n'
#      , xlab = "Brook Char Linkage Group", ylab = "Fst",
#      las = 1, cex.axis = 0.8)
# # Add grey boxes for LGs
# rect(xleft = left.grey[0:21], xright= right.grey[0:21],
#      ybottom = 0, ytop = 11,
#      col = c("lightgrey"),
#      border=NA)
# 
# # Find position for labels
# position <- ((right.grey[1:21] - left.grey[1:21])/2) + left.grey[1:21]
# 
# # Add axis and points
# axis(side = 1, at = position, labels = seq(from = 1, to = 41, by = 2)
#      , cex.axis = 0.8 )
# points(gwas.fst$totpos, gwas.fst$fst, type = "p", cex = 0.8)
# 
# # save as 10 x 6


head(coverage.df)



pdf(file = "barplot_poly_per_region_by_LG.pdf", width = 20, height = 5)
barplot(coverage.df$plotting.poly.sum
        , ylab = "Variants per window (*1000)"
        , las = 1
        )

abline(v = bins.df$end.bin)
axis(side = 1, at = bins.df$label.pos, labels = bins.df$chr)

dev.off()

# Mean by LG
LG_mean.df <- aggregate(coverage.df$poly.sum, by=list(coverage.df$LG), FUN=mean)

LG_sum.df <- aggregate(coverage.df$poly.sum, by=list(coverage.df$LG), FUN=sum)


boxplot(LG_sum.df$x ~ LG_sum.df$Group.1)

pdf(file = "boxplot_poly_per_region_by_LG.pdf", width = 20, height = 5)
boxplot(coverage.df$poly.sum ~ coverage.df$LG)
dev.off()

LG_mean.df
LG_sum.df # note: this should be controlled by the length of the chr, otherwise is not a legitimate comparison

# Stat testing
mod1 <- aov(coverage.df$poly.sum ~ coverage.df$LG)
summary(mod1)
tukey.out <- TukeyHSD(x = mod1)
tukey.df <- as.data.frame(tukey.out$`coverage.df$LG`)
tukey.df[grep(pattern = "LG24", x = rownames(tukey.df)), ]
tukey.df[grep(pattern = "LG09", x = rownames(tukey.df)), ]
tukey.df[grep(pattern = "LG11", x = rownames(tukey.df)), ]
tukey.df[grep(pattern = "LG01", x = rownames(tukey.df)), ]

coverage.df$enriched.chr <- NA
coverage.df$enriched.chr[grep(pattern = "LG09|LG11|LG24", x = coverage.df$LG)] <- "yes"
coverage.df$enriched.chr[grep(pattern = "LG09|LG11|LG24", x = coverage.df$LG, invert = T)] <- "no"

table(coverage.df$enriched.chr)

boxplot(coverage.df$poly.sum ~ coverage.df$enriched.chr)
mod2 <- aov(coverage.df$poly.sum ~ coverage.df$enriched.chr)
summary(mod2)

summary(coverage.df$poly.sum[coverage.df$enriched.chr=="yes"])
summary(coverage.df$poly.sum[coverage.df$enriched.chr=="no"])

table(coverage.df$LG)


# Next steps: 
# normalize to length if doing direct comparisons
# are the window sizes the correct sizes, or should we make them larger?
  # if the goal is to look for 'large swaths', then the window size should probably be larger
# is a Circos plot with histogram for polymorphism level and another track repetitive elements possible? 
   # for this would need to find the repeat content that was used to build the repeat figure

# if make bigger window size
shapiro.test(coverage.df$poly.sum) # data is not normal


#### with advice from deus ex biostats ####
head(coverage.df)
line.of.interest <- coverage.df[which(coverage.df$poly.sum == max(coverage.df$poly.sum)), ]
coverage_except_target.df <- coverage.df[-which(coverage.df$poly.sum == max(coverage.df$poly.sum)), ]

mod <- lm(coverage_except_target.df$poly.sum - line.of.interest$poly.sum ~ 1)
summary(mod)

# Run an iterative analysis for each window in the dataset

slice.of.interest <- NULL; name <- NULL; remainder.df <- NULL; result.list <- list()
mod <- NULL
for(i in 1:nrow(coverage.df)){
  
  # generate a name for the slice
  name <- paste0(coverage.df$LG[i], "_", coverage.df$stop[i])
  
  # extract the slice of interest
  slice.of.interest <- coverage.df[i,]
  
  # Remove the slice of interest from the df
  remainder.df <- coverage.df[-i,]
  
  # Run the model
  mod <- lm(remainder.df$poly.sum - slice.of.interest$poly.sum ~ 1)
  result.list[[name]] <- summary(mod)$coefficients[4]
  
}

result.list
pval_results <- unlist(result.list)

hist(pval_results, breaks = 1000)
hist(coverage.df$poly.sum, breaks = 1000)
boxplot(coverage.df$poly.sum, las = 1, ylab = "per window polymorphism")
points(y = coverage.df$poly.sum[coverage.df$LG=="LG24"])

# Points
stripchart(coverage.df$poly.sum[coverage.df$LG=="LG24"],              # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 4,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

stripchart(coverage.df$poly.sum[coverage.df$LG=="LG09"],              # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 3,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

stripchart(coverage.df$poly.sum[coverage.df$LG=="LG11"],              # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 5,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over



#### OUTLIER APPROACH ####
# How many values are present in the dataset? 
length(coverage.df$poly.sum) #930

# Summary statistics of the polymorphism level per window
summary(coverage.df$poly.sum) # median = 906

# Explore the outliers
boxplot.stats(coverage.df$poly.sum)
# note: confidence interval is 873.4 - 938.6

# Identify the upper outliers: 
boxplot.stats(coverage.df$poly.sum)$out[boxplot.stats(coverage.df$poly.sum)$out > median(coverage.df$poly.sum)]
# note: they are all upper outliers, no lower outliers identified

outlier_polymorphism.vec <- boxplot.stats(coverage.df$poly.sum)$out

min(outlier_polymorphism.vec) # everything over 2341 polymorphism is an outlier
boxplot(coverage.df$poly.sum)
abline(h = min(outlier_polymorphism.vec))

coverage.df[which(coverage.df$poly.sum > min(outlier_polymorphism.vec)), ]

# Inspect this manually to see regions of consecutive elevated polymorphism above the outlier level















