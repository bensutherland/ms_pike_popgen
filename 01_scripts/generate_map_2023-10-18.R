# Plot the sampling locations of the northern pike sampling
# B. Sutherland (2023-10-18)

# Clear space
#rm(list=ls())

# Load packages
#install.packages("vcfR")
#install.packages("rstudioapi")
# install.packages("ggrepel")
# install.packages("maps")
# install.packages("mapdata")
library("ggrepel")
library("maps")
library("mapdata")
library("vcfR")
library("rstudioapi")

# Set working directory to the ms_scallop_popgen repo
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()


#### 01. Provide locations ####
# Manually enter GPS coordinates of each sampling location (long, lat)
CHT <- c(-148.86032, 64.98396)
HOO <- c(-135.13209, 61.51156)
PAL <- c(-133.57592, 59.43708)
CHA <- c(-120.97835, 56.32853)
CAS <- c(-117.65344, 49.31538)
WHI <- c( -95.17243, 49.80051)
SLA <- c( -76.09785, 44.24787)
HCK <- c( -74.83359, 40.84155)


# Make into df
loc_global.df <- rbind(CHT, HOO, PAL, CHA, CAS, WHI, SLA, HCK)
rm(CHT, HOO, PAL, CHA, CAS, WHI, SLA, HCK)
loc_global.df <- as.data.frame(x = loc_global.df, stringsAsFactors = F)
colnames(loc_global.df) <- c("lon", "lat")
loc_global.df

# Extract relevant info
visit.x <- loc_global.df$lon
visit.y <- loc_global.df$lat


#### 2. Mapping ####
# Option 1. maps package
## < probably too detailed for global, use ggplot below >
# pdf(file = "maps_version_global.pdf")
# map("worldHires"
#     #, regions = "."
#     #, xlim=c(-130,-122)
#     #, wrap = c(0,360)
#     , ylim=c(15, 75)
#     , col="gray70"
#     , fill=TRUE
# )  
# 
# points(loc_global.df$lon, loc_global.df$lat, pch=19, col="red", cex=1)  #plot my sample sites
# text(x = loc_global.df$lon, y = loc_global.df$lat, labels =  rownames(loc_global.df), adj = -0.2, cex = 1.2)
# dev.off()


# Option 2. ggplot
# Set variables
ylim.upper <- 60
ylim.lower <- 40
xlim.left <- -180
xlim.right <- -80

map_trimmed <- NULL ; map_plot <- NULL
map_trimmed <- borders(database = "world"
                       , xlim = c(xlim.left, xlim.right)
                       , ylim=c(ylim.lower, ylim.upper)
                       , colour="grey50"
                       , fill="gray70"
)

map_plot <- ggplot() + map_trimmed
map_plot

## Add the NACD
nacd.df <- read.csv(file = "00_archive/NACD_GPS_coordinates_2023-10-18.csv")
head(nacd.df)
tail(nacd.df)

nacd.df <- nacd.df[grep(pattern = "NACD_-02", x = nacd.df$Site, invert = T),]
nacd.df <- nacd.df[grep(pattern = "NACD_-01", x = nacd.df$Site, invert = T),]
nacd.df <- nacd.df[grep(pattern = "NACD_06", x = nacd.df$Site, invert = T),]
nacd.df <- nacd.df[grep(pattern = "NACD_09", x = nacd.df$Site, invert = T),]

# Save point 
#map_plot.bck <- map_plot

# Connect locations with single line through all sites 
map_plot <- map_plot + geom_line(
  aes(x = nacd.df$Longitude
      , y = nacd.df$Latitude
  )
  , col = "grey20"
  , linetype = 2
)
map_plot


# Add sample location points
map_plot <- map_plot + geom_point(aes(x=visit.x, y=visit.y), color="black", size=2) 
map_plot


# Add sample location labels
map_plot <- map_plot + 
  geom_label_repel(
    
    # # Add parameters
    # direction = "x"
    # , vjust = 1
    # , hjust = 1
    
    # Add aesthetics
    aes(x = loc_global.df$lon
          , y = loc_global.df$lat
          , label = rownames(loc_global.df)
    )
    
  )

map_plot

# Add NACD start label
map_plot <- map_plot + 
  #geom_label_repel(
  geom_text( 
    aes(  
        #x = nacd.df$Longitude[nrow(nacd.df)-1]
          x = nacd.df$Longitude[1]
        #, y =  nacd.df$Latitude[nrow(nacd.df)-1]
        , y =  nacd.df$Latitude[1]
        , label = "NACD"
        #, label.size = 0.5
        #, colour="#06414b"
        #, segment.colour="black"
        
        
        )
    
  )

map_plot

# Save out
pdf(file = "03_results/northern_pike_map_figure.pdf", width = 9, height = 5
)
map_plot
dev.off()

