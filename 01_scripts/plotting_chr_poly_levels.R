# Visualizing the differences in chr polymorphism levels across and among chromosomes
# Ben J. G. Sutherland, 2023-11-01

# Clear space
#rm(list=ls())

# Load packages
#install.packages("ggplot2")
library("ggplot2")

# Set working directory to the ms_scallop_popgen repo
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# Set variables
options(scipen = 9999999)
coverage.FN <- "02_input_data/coverage.txt"

#### 00. Read in data ####
coverage.df <- read.delim(file = coverage.FN, header = F)
head(coverage.df, n = 10)
colnames(coverage.df) <- c("LG", "start", "stop", "poly.sum")
dim(coverage.df)
str(coverage.df)
coverage.df$start <- as.numeric(coverage.df$start)
coverage.df$stop  <- as.numeric(coverage.df$stop)

# Sort by chr
coverage.df <- coverage.df[with(coverage.df, order(coverage.df$LG, coverage.df$start)), ]
unique(coverage.df$LG)
head(coverage.df)


#### 01. Polymorphism by chr length ####
# Determine the length of each chr
len.chr = NULL; len.chr.list <- list()
for(i in 1:length(unique(coverage.df$LG))){

  chr.o.i <- unique(coverage.df$LG)[i]
  print(chr.o.i)

  len.chr.list[[chr.o.i]] <- max(coverage.df[coverage.df$LG==chr.o.i, "stop"])

}

len.chr.df <- as.data.frame(unlist(len.chr.list))
len.chr.df$chr <- rownames(len.chr.df)
colnames(len.chr.df) <- c("len", "chr")
len.chr.df <- len.chr.df[,c("chr", "len")]
head(len.chr.df)
tail(len.chr.df)

# Calculate number of variants per chr
poly_by_chr.df <- aggregate(coverage.df$poly.sum, by=list(LG=coverage.df$LG), FUN=sum)
str(poly_by_chr.df)
head(poly_by_chr.df)
# test
sum(coverage.df[coverage.df$LG=="LG01", "poly.sum"]) # to confirm matches above

# Combine variant counts with chr
head(poly_by_chr.df)
head(len.chr.df)
chr_length_and_cumul_poly.df <- merge(x = len.chr.df, y = poly_by_chr.df, by.x = "chr", by.y = "LG")
head(chr_length_and_cumul_poly.df)
colnames(chr_length_and_cumul_poly.df) <- c("chr", "length", "poly.sum")
head(chr_length_and_cumul_poly.df)

# Convert bp to kbp
chr_length_and_cumul_poly.df$length_per.kb <- chr_length_and_cumul_poly.df$length / 1000
head(chr_length_and_cumul_poly.df) # This will be used below for plotting


#### 02. Plot windows (binned variants) across chr ####
# How many bins of variants per chr? 
bins <- table(coverage.df$LG)
bins.df <- as.data.frame(names(bins))
colnames(bins.df) <- "chr"
bins.df$bins <- as.numeric(bins)
bins.df$cumsum <- cumsum(bins.df$bins) # cumulate total number of bins

head(bins.df)

# Determine bin numbers where chr start, stop, or midway for labels
bins.df$start.bin <- NA; bins.df$end.bin <- NA; bins.df$label.pos <- NA
for(i in 1:nrow(bins.df)){
  
  if(i == 1){
    
    # The starting chr starts at zero
    bins.df$start.bin[i] <- 0
    
  }else{
    
    # Any subsequent chr starts at the end of the previous cumulative bin count
    bins.df$start.bin[i] <- bins.df$cumsum[i-1]
    
  }
  
  # End of bins is at the end of the cumulative sum per chr
  bins.df$end.bin[i] <- bins.df$cumsum[i]
  
  # Labels go in the middle of each start and stop per chr
  bins.df$label.pos[i] <- ((bins.df$end.bin[i] - bins.df$start.bin[i])/2) + bins.df$start.bin[i]

}

head(bins.df)

# Create alternating colour vector for plotting
bins.df$plot.colour <- "black"
bins.df$plot.colour[c(FALSE, TRUE)] <- "darkgrey"
head(bins.df)

# Inspect df to be used
head(coverage.df)
head(bins.df)

# Bring colour from the bins df into the coverage df
coverage.df$plot.color <- NA
for(i in 1:length(coverage.df$plot.color)){
  
  coverage.df$plot.color[i] <- bins.df[bins.df$chr %in% coverage.df$LG[i], "plot.colour"]
  
}

head(coverage.df, n = 10)


#### 03. Identify bins that have outlier variant counts ####
# How many bins were generated? 
length(coverage.df$poly.sum) #930

# Summary statistics of the polymorphism level per window
summary(coverage.df$poly.sum) # median = 906

# Explore the outliers
boxplot.stats(coverage.df$poly.sum)
# note: confidence interval is 873.4 - 938.6

# Identify the upper outliers: 
boxplot.stats(coverage.df$poly.sum)$out[boxplot.stats(coverage.df$poly.sum)$out > median(coverage.df$poly.sum)]
#  note: they are all upper outliers, no lower outliers identified

outlier_polymorphism.vec <- boxplot.stats(coverage.df$poly.sum)$out

min(outlier_polymorphism.vec) # everything over 2341 polymorphism is an outlier
outlier_cutoff.val <- min(outlier_polymorphism.vec) # everything over 2341 polymorphism is an outlier

boxplot(coverage.df$poly.sum, las = 1)
abline(h = outlier_cutoff.val, lty = 2)

# What bins are above the outlier level? 
head(coverage.df[which(coverage.df$poly.sum > min(outlier_polymorphism.vec)), ])
outlier_segments.df <- coverage.df[which(coverage.df$poly.sum > min(outlier_polymorphism.vec)), ]

# How many from each chr? 
table(coverage.df[which(coverage.df$poly.sum > min(outlier_polymorphism.vec)), "LG"])

# Export results as text files
write.table(x = outlier_segments.df, file = "03_results/elevated_variant_windows.txt", quote = F
            , sep = "\t", row.names = F, col.names = T
            ) # outliers only

write.table(x = coverage.df, file = "03_results/all_variant_windows.txt", quote = F
            , sep = "\t", row.names = F, col.names = T
            ) # all data

# Note: inspect this manually to see regions of consecutive elevated polymorphism above the outlier level
head(coverage.df)

coverage.df <- coverage.df[with(coverage.df, order(coverage.df$LG, coverage.df$start)), ]
head(coverage.df)
tail(coverage.df)
coverage.df$plot.order <- seq(1:nrow(coverage.df))

#### 04. Plotting ####
pdf(file = "03_results/variants_per_chr.pdf", width = 11, height = 6)
par(mfrow=c(2,1), mar = c(4.5,4.1,3, 2.1))
# Plot variant sum per bin across the chr
#pdf(file = "03_results/plot_SNP_sum_per_window.pdf", width = 15, height = 4)
plot(x = coverage.df$plot.order, y = coverage.df$poly.sum
     , xaxt = "n"
     , xlab = "Chromosome"
     , ylab = "variants per window"
     , las = 1
     , pch = 16
     , col = coverage.df$plot.color
     , cex = 0.8
     )
abline(v = bins.df$end.bin, lty = 2)
abline(v = bins.df$start.bin, lty = 2)
axis(side = 1, at = bins.df$label.pos, labels = gsub(pattern = "LG", replacement = "", x = bins.df$chr)
     , tick = T
     , las = 2
     )

abline(h = outlier_cutoff.val, lty = 3)

#dev.off()


# Plot mean number of variants per kbp per chr
#pdf(file = "03_results/plot_mean_variant_count_per_window_per_chr.pdf", width = 15, height = 4)
plot(chr_length_and_cumul_poly.df$poly.sum / chr_length_and_cumul_poly.df$length_per.kb
     , xaxt = "n"
     , ylab = "avg. variant per kb"
     , las = 1
     , xlab = "Chromosome"
     , pch = 16
     , ylim = c(0,3)
)
axis(side = 1, at = seq(1:nrow(chr_length_and_cumul_poly.df))
     , labels = gsub(pattern = "LG", replacement = "", x = chr_length_and_cumul_poly.df$chr)
     , las = 2)

#dev.off()
dev.off()







#### OLD CODE ####

# pdf(file = "barplot_poly_per_region_by_LG.pdf", width = 20, height = 5)
# barplot(coverage.df$plotting.poly.sum
#         , ylab = "Variants per window (*1000)"
#         , las = 1
#         )
# 
# abline(v = bins.df$end.bin)
# axis(side = 1, at = bins.df$label.pos, labels = bins.df$chr)
# 
# dev.off()

# # Mean variants per bin per chr (OLD)
# LG_mean.df <- aggregate(coverage.df$poly.sum, by=list(coverage.df$LG), FUN=mean)
# LG_sum.df <- aggregate(coverage.df$poly.sum, by=list(coverage.df$LG), FUN=sum)
# #boxplot(LG_sum.df$x ~ LG_sum.df$Group.1)

# # Boxplot of distribution of variants per bin per chr (OLD)
# pdf(file = "boxplot_poly_per_region_by_LG.pdf", width = 20, height = 5)
# boxplot(coverage.df$poly.sum ~ coverage.df$LG)
# dev.off()

# LG_mean.df
# LG_sum.df # note: this should be controlled by the length of the chr, otherwise is not a legitimate comparison

# #### Stat testing of num variants per chr (still needs to be normalized by length) ###
# mod1 <- aov(coverage.df$poly.sum ~ coverage.df$LG)
# summary(mod1)
# tukey.out <- TukeyHSD(x = mod1)
# tukey.df <- as.data.frame(tukey.out$`coverage.df$LG`)
# tukey.df[grep(pattern = "LG24", x = rownames(tukey.df)), ]
# tukey.df[grep(pattern = "LG09", x = rownames(tukey.df)), ]
# tukey.df[grep(pattern = "LG11", x = rownames(tukey.df)), ]
# tukey.df[grep(pattern = "LG01", x = rownames(tukey.df)), ]
# 
# coverage.df$enriched.chr <- NA
# coverage.df$enriched.chr[grep(pattern = "LG09|LG11|LG24", x = coverage.df$LG)] <- "yes"
# coverage.df$enriched.chr[grep(pattern = "LG09|LG11|LG24", x = coverage.df$LG, invert = T)] <- "no"
# 
# table(coverage.df$enriched.chr)
# 
# boxplot(coverage.df$poly.sum ~ coverage.df$enriched.chr)
# mod2 <- aov(coverage.df$poly.sum ~ coverage.df$enriched.chr)
# summary(mod2)
# 
# summary(coverage.df$poly.sum[coverage.df$enriched.chr=="yes"])
# summary(coverage.df$poly.sum[coverage.df$enriched.chr=="no"])
# 
# table(coverage.df$LG)
# #### /END/ Stat testing of num variants per chr (still needs to be normalized by length) ###

# if make bigger window size
#shapiro.test(coverage.df$poly.sum) # data is not normal


# #### with advice from deus ex biostats ####
# head(coverage.df)
# line.of.interest <- coverage.df[which(coverage.df$poly.sum == max(coverage.df$poly.sum)), ]
# coverage_except_target.df <- coverage.df[-which(coverage.df$poly.sum == max(coverage.df$poly.sum)), ]
# 
# mod <- lm(coverage_except_target.df$poly.sum - line.of.interest$poly.sum ~ 1)
# summary(mod)
# 
# # Run an iterative analysis for each window in the dataset
# 
# slice.of.interest <- NULL; name <- NULL; remainder.df <- NULL; result.list <- list()
# mod <- NULL
# for(i in 1:nrow(coverage.df)){
#   
#   # generate a name for the slice
#   name <- paste0(coverage.df$LG[i], "_", coverage.df$stop[i])
#   
#   # extract the slice of interest
#   slice.of.interest <- coverage.df[i,]
#   
#   # Remove the slice of interest from the df
#   remainder.df <- coverage.df[-i,]
#   
#   # Run the model
#   mod <- lm(remainder.df$poly.sum - slice.of.interest$poly.sum ~ 1)
#   result.list[[name]] <- summary(mod)$coefficients[4]
#   
# }
# 
# result.list
# pval_results <- unlist(result.list)
# 
# hist(pval_results, breaks = 1000)
# hist(coverage.df$poly.sum, breaks = 1000)
# boxplot(coverage.df$poly.sum, las = 1, ylab = "per window polymorphism")
# #points(y = coverage.df$poly.sum[coverage.df$LG=="LG24"])
# 
# # Points
# stripchart(coverage.df$poly.sum[coverage.df$LG=="LG24"],              # Data
#            method = "jitter", # Random noise
#            pch = 19,          # Pch symbols
#            col = 4,           # Color of the symbol
#            vertical = TRUE,   # Vertical mode
#            add = TRUE)        # Add it over
# 
# stripchart(coverage.df$poly.sum[coverage.df$LG=="LG09"],              # Data
#            method = "jitter", # Random noise
#            pch = 19,          # Pch symbols
#            col = 3,           # Color of the symbol
#            vertical = TRUE,   # Vertical mode
#            add = TRUE)        # Add it over
# 
# stripchart(coverage.df$poly.sum[coverage.df$LG=="LG11"],              # Data
#            method = "jitter", # Random noise
#            pch = 19,          # Pch symbols
#            col = 5,           # Color of the symbol
#            vertical = TRUE,   # Vertical mode
#            add = TRUE)        # Add it over

