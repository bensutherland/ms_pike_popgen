# Generate multi-panel plot for publication
# B. Sutherland (2023-11-17)

# Clear space
#rm(list=ls())

# Load packages
#install.packages("rstudioapi")
#install.packages("ggpubr")
#install.packages("ggside")
# utils::install.packages("ggside")
# utils::install.packages("ggstatsplot")
#install.packages("ggstatsplot")
library("rstudioapi")
require("ggpubr")
library("tidyr")
library("vcfR")
library("rstudioapi")
library("dplyr")
library("ggplot2")
library("ggside")
library("ggstatsplot")

# Set working directory to the ms_scallop_popgen repo
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# load(file = "03_results/sample_map.Rdata")
# load(file = "03_results/NACD_df.Rdata")
load(file = "03_results/mapping_data_with_MIN.Rdata")

load(file = "03_results/per_ind_genos.Rdata")

figure <- ggarrange(map_plot, variants_per_indiv.plot
                    , labels = c("A", "B")
                    , ncol = 1, nrow = 2
)


pdf(file = "03_results/figure_map_and_geno_counts.pdf", width = 7.4, height = 7)
print(figure)
dev.off()
