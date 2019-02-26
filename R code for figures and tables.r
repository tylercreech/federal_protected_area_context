# GRAPHICAL ANALYSIS OF BIOLOGICAL VARIABLES

library(raster)
library(sp)
library(sf)
library(rgeos)
library(rgdal)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(rvest)
library(maptools)
library(ggplot2)
library(RColorBrewer)
library(rasterVis)
library(gridExtra)
library(diveRsity)
library(grid)
library(gridExtra)
library(legendMap)
library(fasterize)


####################################################################################################
#### DATA PREP
####################################################################################################

infolder <- "H:/My Drive/MonumentData/Generated Data"  # set folder holding input data

# load RData file containing zonal stats in dataframe called PA_zonal.df
PA_zonal.df <- read.csv(paste0(infolder, "/PA_zonal_stats_post1996_lower48_14Feb2019.csv"), head=TRUE, stringsAsFactors = FALSE)


# load spatial data
PA <- st_read(paste(infolder, "/post1996_PAs_lower48_14Feb2019.shp", sep=""))  # use this version of the PA shapefile that has duplicate Unit Names corrected - maps will be incorrect otherwise
# rename PA attributes (st_write abbreviates these to match ESRI character limitations)
names(PA) <- c("UnitName", "InReview", "CurDesType","CurDesAuth", "OriDesAuth", "CurDesYear", "AntiqYear", "area_m2", "area_ac", "DesMode", "clipped_area_m2", "clipped_area_ac", "clipped_fraction", "geometry")
rich.bird <- raster(paste(infolder, "/rich.bird.tif", sep=""))
rich.mammal <- raster(paste(infolder, "/rich.mammal.tif", sep=""))
rich.tree <- raster(paste(infolder, "/rich.tree.tif", sep=""))
rich.reptile <- raster(paste(infolder, "/rich.reptile.tif", sep=""))
rich.natserv <- raster(paste(infolder, "/natserv.tif", sep=""))
rich.fish <- st_read(paste(infolder, "/rich.fish.shp", sep=""))
rich.amphib <- st_read(paste(infolder, "/rich.amphib.shp", sep=""))
natlandcover <- raster(paste(infolder, "/natlandcover.tif", sep=""))
climate <- raster(paste(infolder, "/climate_clipped.tif", sep=""))
fedlands <- st_read(paste(infolder, "/fedlandscrop.shp", sep=""))
states <- st_read(paste(infolder, "/states2.shp", sep=""))


### MANDATORY FILTERING

# get rid of PAs that are too impervious or are redesignated PAs
discards <- PA_zonal.df$UnitName[PA_zonal.df$keep=="no"]
PA <- PA[-which(PA$UnitName %in% discards),]
PA_zonal.df <- PA_zonal.df[which(PA_zonal.df$keep=="yes"),]



### OPTIONAL FILTERING

# subset zonal data to eastern US
westernstates <- c("Washington", "Oregon", "California", "Nevada", "Idaho", "Montana", "Wyoming", "Utah", "Arizona", "Colorado", "New Mexico")
PA_zonal.df <- PA_zonal.df[-which(PA_zonal.df$state.majority %in% westernstates),]

# subset zonal data to western US
westernstates <- c("Washington", "Oregon", "California", "Nevada", "Idaho", "Montana", "Wyoming", "Utah", "Arizona", "Colorado", "New Mexico")
PA_zonal.df <- PA_zonal.df[which(PA_zonal.df$state.majority %in% westernstates),]




###################################################################################################
### MAP OF PAs BY DESIGNATION MODE WITH SIZE HISTOGRAM AS INSERT
###################################################################################################


### Map portion showing PAs by designation type
fedlands.raster <- fasterize(sf=fedlands, raster=climate, field=NULL, fun="last")  # rasterize fedlands layer to allow efficient plotting
states.wgs84 <- st_transform(states, "+init=epsg:4326")  # transform CRS to WGS84 (necessary to use legendMap for scale bar)
PA.wgs84 <- st_transform(PA, "+init=epsg:4326")
fedlands.raster.wgs84 <- projectRaster(fedlands.raster, crs="+init=epsg:4326", method="ngb")
fedlands.df <- as.data.frame(as(fedlands.raster.wgs84, "SpatialPixelsDataFrame"))  # convert to dataframe to allow plotting using geom_tile
colnames(fedlands.df) <- c("value", "x", "y")

PAmap <- ggplot() +
  geom_sf(data=states.wgs84, fill="grey90", color=NA) +  # light grey background for lower 48 states
  geom_tile(data=fedlands.df, aes(x=x, y=y, fill="zero")) +  # federal lands in darker grey
  geom_sf(data=states.wgs84, fill=NA, color="white") +  # state outlines in darkest grey
  geom_sf(data=PA.wgs84[which(PA.wgs84$DesMode=="Congress"),], aes(fill="first"), color=NA, alpha=0.5) +  # protected areas on top
  geom_sf(data=PA.wgs84[which(PA.wgs84$DesMode=="President"),], aes(fill="second"), color=NA, alpha=0.5) +  # protected areas on top
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        panel.grid.major = element_line(colour = "white")) +
  scale_fill_manual(name="Designation\nmode", values=c("first"="red","second"="blue",
        "zero"="gray70"), labels=c("CPA","PPA","Federal land")) +
  theme(legend.justification=c(1,0), legend.position="bottom") +   # put legend in top right corner
  scale_bar(lon=-75, lat=22, distance_lon=500, distance_lat=50, 
            distance_legend=150, dist_unit="km", rec_fill="white", 
            rec_colour="black", rec2_fill="black", rec2_colour="black", 
            legend_colour="black", legend_size=3, orientation=TRUE, 
            arrow_length=500, arrow_distance=300, arrow_north_size=4)

# Density plot of PA size distribution
bandwidth.multiplier <- 1   # adjust bandwidth (how much smoothing of histogram occurs - 1 is default)
pointsize <- 2   # adust size of points in error bars
heightfactor  <- 50  # adjust height of tick marks on end of error bars  (height = yrange of plot / heightfactor )
PA_zonal.df$clipped_area_km2 <- PA_zonal.df$clipped_area_m2 / 1e+06
PA_zonal.df$log_clipped_area_km2 <- log(PA_zonal.df$clipped_area_km2)
attach(PA_zonal.df)
size.mean.ppa <- log(mean(clipped_area_km2[which(DesMode=="President")], na.rm=TRUE))
size.mean.cpa <- log(mean(clipped_area_km2[which(DesMode=="Congress")], na.rm=TRUE))
size.sd.ppa <- log(sd(clipped_area_km2[which(DesMode=="President")], na.rm=TRUE))
size.sd.cpa <- log(sd(clipped_area_km2[which(DesMode=="Congress")], na.rm=TRUE))
size <- ggplot() +
   geom_density(data=PA_zonal.df, aes(x=log_clipped_area_km2, fill=DesMode, color=DesMode), alpha=0.35, size=1) + 
   labs(x="Log of PA area (sq. km)", y="Density") +
   scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
   scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
   theme(legend.position="none",
                   rect = element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   panel.background = element_blank(),
                   axis.line.x = element_line(),
                   axis.line.y = element_line())
size.yrange <- layer_scales(size)$y$range$range[2]  # get y range of density plot
size.errorbar <- size +   # add error bars to plot
   geom_errorbarh(aes(y=-size.yrange/10, xmax=size.mean.ppa + size.sd.ppa, xmin=size.mean.ppa - size.sd.ppa, height=size.yrange/heightfactor), color="blue", size=1) +
   geom_point(aes(y=-size.yrange/10, x=size.mean.ppa), color="blue", size=pointsize) +
   geom_errorbarh(aes(y=-size.yrange/15, xmax=size.mean.cpa + size.sd.cpa, xmin=size.mean.cpa - size.sd.cpa, height=size.yrange/heightfactor), color="red", size=1) +
   geom_point(aes(y=-size.yrange/15, x=size.mean.cpa), color="red", size=pointsize)


# Put map and density plot in same window
grid.newpage()
large <- viewport(width = 1, height = 1, x = 0.5, y = 0.5)  # the larger map
small <- viewport(width = 0.4, height = 0.4, x = 0.22, y = 0.21)  # the inset in upper right
print(PAmap, vp = large)
print(size.errorbar, vp = small)






###################################################################################################
### DENSITY PLOTS OF ECOLOGICAL VARIABLES BY DESIGNATION MODE
###################################################################################################


# ADD bandwidth adjustment here
bandwidth.multiplier <- 1   # adjust bandwidth (how much smoothing of histogram occurs - 1 is default)
pointsize <- 2   # adust size of points in error bars
heightfactor  <- 50  # adjust height of tick marks on end of error bars  (height = yrange of plot / heightfactor )
stat <- "mean"   # choose within-PA summary stat; "mean", "min", or "max"
statname <- "Mean"  # write out statistic name as you want it to appear in plot legend (full name, capitalized first letter:  "Mean", "Minimum", or "Maximum)
attach(PA_zonal.df)

# bird richness
bird.vals <- get(paste0(stat,".rich.bird"))  # get mean and sd of value for each mode
bird.mean.ppa <- mean(bird.vals[which(DesMode=="President")], na.rm=TRUE)
bird.mean.cpa <- mean(bird.vals[which(DesMode=="Congress")], na.rm=TRUE)
bird.sd.ppa <- sd(bird.vals[which(DesMode=="President")], na.rm=TRUE)
bird.sd.cpa <- sd(bird.vals[which(DesMode=="Congress")], na.rm=TRUE)
r1 <- ggplot() +   # make density plot
  geom_density(data=PA_zonal.df, aes(x=get(paste0(stat,".rich.bird")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x=paste0(statname, " # species"), y="Density") +
  ggtitle("Bird richness") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
bird.yrange <- layer_scales(r1)$y$range$range[2]  # get y range of density plot
r1.errorbar <- r1 +   # add error bars to plot
  geom_errorbarh(aes(y=-bird.yrange/10, xmax=bird.mean.ppa + bird.sd.ppa, xmin=bird.mean.ppa - bird.sd.ppa, height=bird.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-bird.yrange/10, x=bird.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-bird.yrange/13, xmax=bird.mean.cpa + bird.sd.ppa, xmin=bird.mean.cpa - bird.sd.cpa, height=bird.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-bird.yrange/13, x=bird.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r1)$y$range$range[2]/20,layer_scales(r1)$y$range$range[2]/20))

# mammal richness
mammal.vals <- get(paste0(stat,".rich.mammal"))  # get mean and sd of value for each mode
mammal.mean.ppa <- mean(mammal.vals[which(DesMode=="President")], na.rm=TRUE)
mammal.mean.cpa <- mean(mammal.vals[which(DesMode=="Congress")], na.rm=TRUE)
mammal.sd.ppa <- sd(mammal.vals[which(DesMode=="President")], na.rm=TRUE)
mammal.sd.cpa <- sd(mammal.vals[which(DesMode=="Congress")], na.rm=TRUE)
r2 <- ggplot() +    
  geom_density(data=PA_zonal.df, aes(x=get(paste0(stat,".rich.mammal")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x=paste0(statname, " # species"), y="Density") +
  ggtitle("Mammal richness") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
mammal.yrange <- layer_scales(r2)$y$range$range[2]  # get y range of density plot
r2.errorbar <- r2 +   # add error bars to plot
  geom_errorbarh(aes(y=-mammal.yrange/10, xmax=mammal.mean.ppa + mammal.sd.ppa, xmin=mammal.mean.ppa - mammal.sd.ppa, height=mammal.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-mammal.yrange/10, x=mammal.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-mammal.yrange/13, xmax=mammal.mean.cpa + mammal.sd.ppa, xmin=mammal.mean.cpa - mammal.sd.cpa, height=mammal.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-mammal.yrange/13, x=mammal.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r2)$y$range$range[2]/20,layer_scales(r2)$y$range$range[2]/20))

# fish richness
fish.vals <- get(paste0(stat,".rich.fish"))  # get mean and sd of value for each mode
fish.mean.ppa <- mean(fish.vals[which(DesMode=="President")], na.rm=TRUE)
fish.mean.cpa <- mean(fish.vals[which(DesMode=="Congress")], na.rm=TRUE)
fish.sd.ppa <- sd(fish.vals[which(DesMode=="President")], na.rm=TRUE)
fish.sd.cpa <- sd(fish.vals[which(DesMode=="Congress")], na.rm=TRUE)
r3 <- ggplot() +    
  geom_density(data=PA_zonal.df, aes(x=get(paste0(stat,".rich.fish")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x=paste0(statname, " # species"), y="Density") +
  ggtitle("Fish richness") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  theme(legend.justification=c(1,1), legend.position=c(1,1)) + # put legend in top right corner
  theme(legend.text=element_text(size=8), legend.title=element_text(size=8))
fish.yrange <- layer_scales(r3)$y$range$range[2]  # get y range of density plot
r3.errorbar <- r3 +   # add error bars to plot
  geom_errorbarh(aes(y=-fish.yrange/10, xmax=fish.mean.ppa + fish.sd.ppa, xmin=fish.mean.ppa - fish.sd.ppa, height=fish.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-fish.yrange/10, x=fish.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-fish.yrange/13, xmax=fish.mean.cpa + fish.sd.ppa, xmin=fish.mean.cpa - fish.sd.cpa, height=fish.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-fish.yrange/13, x=fish.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r3)$y$range$range[2]/20,layer_scales(r3)$y$range$range[2]/20))

# amphibian richness
amphib.vals <- get(paste0(stat,".rich.amphib"))  # get mean and sd of value for each mode
amphib.mean.ppa <- mean(amphib.vals[which(DesMode=="President")], na.rm=TRUE)
amphib.mean.cpa <- mean(amphib.vals[which(DesMode=="Congress")], na.rm=TRUE)
amphib.sd.ppa <- sd(amphib.vals[which(DesMode=="President")], na.rm=TRUE)
amphib.sd.cpa <- sd(amphib.vals[which(DesMode=="Congress")], na.rm=TRUE)
r4 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0(stat,".rich.amphib")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x=paste0(statname, " # species"), y="Density") +
  ggtitle("Amphibian richness") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
amphib.yrange <- layer_scales(r4)$y$range$range[2]  # get y range of density plot
r4.errorbar <- r4 +   # add error bars to plot
  geom_errorbarh(aes(y=-amphib.yrange/10, xmax=amphib.mean.ppa + amphib.sd.ppa, xmin=amphib.mean.ppa - amphib.sd.ppa, height=amphib.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-amphib.yrange/10, x=amphib.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-amphib.yrange/13, xmax=amphib.mean.cpa + amphib.sd.ppa, xmin=amphib.mean.cpa - amphib.sd.cpa, height=amphib.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-amphib.yrange/13, x=amphib.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r4)$y$range$range[2]/20,layer_scales(r4)$y$range$range[2]/20))

# reptile richness
reptile.vals <- get(paste0(stat,".rich.reptile"))  # get mean and sd of value for each mode
reptile.mean.ppa <- mean(reptile.vals[which(DesMode=="President")], na.rm=TRUE)
reptile.mean.cpa <- mean(reptile.vals[which(DesMode=="Congress")], na.rm=TRUE)
reptile.sd.ppa <- sd(reptile.vals[which(DesMode=="President")], na.rm=TRUE)
reptile.sd.cpa <- sd(reptile.vals[which(DesMode=="Congress")], na.rm=TRUE)
r5 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0(stat,".rich.reptile")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x=paste0(statname, " # species"), y="Density") +
  ggtitle("Reptile richness") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
reptile.yrange <- layer_scales(r5)$y$range$range[2]  # get y range of density plot
r5.errorbar <- r5 +   # add error bars to plot
  geom_errorbarh(aes(y=-reptile.yrange/10, xmax=reptile.mean.ppa + reptile.sd.ppa, xmin=reptile.mean.ppa - reptile.sd.ppa, height=reptile.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-reptile.yrange/10, x=reptile.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-reptile.yrange/13, xmax=reptile.mean.cpa + reptile.sd.ppa, xmin=reptile.mean.cpa - reptile.sd.cpa, height=reptile.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-reptile.yrange/13, x=reptile.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r5)$y$range$range[2]/20,layer_scales(r5)$y$range$range[2]/20))

# tree richness
tree.vals <- get(paste0(stat,".rich.tree"))  # get mean and sd of value for each mode
tree.mean.ppa <- mean(tree.vals[which(DesMode=="President")], na.rm=TRUE)
tree.mean.cpa <- mean(tree.vals[which(DesMode=="Congress")], na.rm=TRUE)
tree.sd.ppa <- sd(tree.vals[which(DesMode=="President")], na.rm=TRUE)
tree.sd.cpa <- sd(tree.vals[which(DesMode=="Congress")], na.rm=TRUE)
r6 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0(stat,".rich.tree")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x=paste0(statname, " # species"), y="Density") +
  ggtitle("Tree richness") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
tree.yrange <- layer_scales(r6)$y$range$range[2]  # get y range of density plot
r6.errorbar <- r6 +   # add error bars to plot
  geom_errorbarh(aes(y=-tree.yrange/10, xmax=tree.mean.ppa + tree.sd.ppa, xmin=tree.mean.ppa - tree.sd.ppa, height=tree.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-tree.yrange/10, x=tree.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-tree.yrange/13, xmax=tree.mean.cpa + tree.sd.ppa, xmin=tree.mean.cpa - tree.sd.cpa, height=tree.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-tree.yrange/13, x=tree.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r6)$y$range$range[2]/20,layer_scales(r6)$y$range$range[2]/20))

# G1/G2 richness
natserv.vals <- get(paste0(stat,".rich.natserv"))  # get mean and sd of value for each mode
natserv.mean.ppa <- mean(natserv.vals[which(DesMode=="President")], na.rm=TRUE)
natserv.mean.cpa <- mean(natserv.vals[which(DesMode=="Congress")], na.rm=TRUE)
natserv.sd.ppa <- sd(natserv.vals[which(DesMode=="President")], na.rm=TRUE)
natserv.sd.cpa <- sd(natserv.vals[which(DesMode=="Congress")], na.rm=TRUE)
r7 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0(stat,".rich.natserv")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x=paste0(statname, " rarity-weighted richness"), y="Density") +
  ggtitle("G1 & G2 species richness") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
natserv.yrange <- layer_scales(r7)$y$range$range[2]  # get y range of density plot
r7.errorbar <- r7 +   # add error bars to plot
  geom_errorbarh(aes(y=-natserv.yrange/10, xmax=natserv.mean.ppa + natserv.sd.ppa, xmin=natserv.mean.ppa - natserv.sd.ppa, height=natserv.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-natserv.yrange/10, x=natserv.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-natserv.yrange/13, xmax=natserv.mean.cpa + natserv.sd.ppa, xmin=natserv.mean.cpa - natserv.sd.cpa, height=natserv.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-natserv.yrange/13, x=natserv.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r7)$y$range$range[2]/80,layer_scales(r7)$y$range$range[2]/80))

# ecological system richness
system.vals <- system.richness.rare  # get mean and sd of value for each mode
system.mean.ppa <- mean(system.vals[which(DesMode=="President")], na.rm=TRUE)
system.mean.cpa <- mean(system.vals[which(DesMode=="Congress")], na.rm=TRUE)
system.sd.ppa <- sd(system.vals[which(DesMode=="President")], na.rm=TRUE)
system.sd.cpa <- sd(system.vals[which(DesMode=="Congress")], na.rm=TRUE)
r8 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(system.richness.rare, fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x="Rarefied richness", y="Density") +
  ggtitle("Ecological system richness") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
system.yrange <- layer_scales(r8)$y$range$range[2]  # get y range of density plot
r8.errorbar <- r8 +   # add error bars to plot
  geom_errorbarh(aes(y=-system.yrange/10, xmax=system.mean.ppa + system.sd.ppa, xmin=system.mean.ppa - system.sd.ppa, height=system.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-system.yrange/10, x=system.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-system.yrange/13, xmax=system.mean.cpa + system.sd.ppa, xmin=system.mean.cpa - system.sd.cpa, height=system.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-system.yrange/13, x=system.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r8)$y$range$range[2]/20,layer_scales(r8)$y$range$range[2]/20))

# climate refugial potential
climate.vals <- get(paste0(stat,".climate"))  # get mean and sd of value for each mode
climate.mean.ppa <- mean(climate.vals[which(DesMode=="President")], na.rm=TRUE)
climate.mean.cpa <- mean(climate.vals[which(DesMode=="Congress")], na.rm=TRUE)
climate.sd.ppa <- sd(climate.vals[which(DesMode=="President")], na.rm=TRUE)
climate.sd.cpa <- sd(climate.vals[which(DesMode=="Congress")], na.rm=TRUE)
r9 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0(stat,".climate")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x=paste0(statname, " refugial index"), y="Density") +
  ggtitle("Climate refugial potential") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
climate.yrange <- layer_scales(r9)$y$range$range[2]  # get y range of density plot
r9.errorbar <- r9 +   # add error bars to plot
  geom_errorbarh(aes(y=-climate.yrange/10, xmax=climate.mean.ppa + climate.sd.ppa, xmin=climate.mean.ppa - climate.sd.ppa, height=climate.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-climate.yrange/10, x=climate.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-climate.yrange/13, xmax=climate.mean.cpa + climate.sd.ppa, xmin=climate.mean.cpa - climate.sd.cpa, height=climate.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-climate.yrange/13, x=climate.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r9)$y$range$range[2]/40,layer_scales(r9)$y$range$range[2]/40))

detach(PA_zonal.df)
#multiplot(r1,r2,r3,r4,r5,r6,r7,r8,r9, cols=3)  # density plots only 
multiplot(r1.errorbar, r2.errorbar, r3.errorbar,r4.errorbar,r5.errorbar,r6.errorbar,r7.errorbar,r8.errorbar,r9.errorbar, cols=3)   # density plots with means and SDs








###################################################################################################
### DENSITY PLOTS OF SOCIOPOLITICAL VARIABLES BY DESIGNATION MODE
###################################################################################################


# MEANS
bandwidth.multiplier <- 1   # adjust bandwidth (how much smoothing of histogram occurs - 1 is default)
pointsize <- 2   # adust size of points in error bars
heightfactor  <- 50  # adjust height of tick marks on end of error bars  (height = yrange of plot / heightfactor )

attach(PA_zonal.df)
bufdist1 <- "10"   # select buffer distance; must be character format, one of these choices: "10","20","50","100","250")
# LCV score
lcv.vals <- get(paste0("mean.LCV.",bufdist1,"kmBuffer"))  # get mean and sd of value for each mode
lcv.mean.ppa <- mean(lcv.vals[which(DesMode=="President")], na.rm=TRUE)
lcv.mean.cpa <- mean(lcv.vals[which(DesMode=="Congress")], na.rm=TRUE)
lcv.sd.ppa <- sd(lcv.vals[which(DesMode=="President")], na.rm=TRUE)
lcv.sd.cpa <- sd(lcv.vals[which(DesMode=="Congress")], na.rm=TRUE)
s1 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0("mean.LCV.",bufdist1,"kmBuffer")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x="Mean LCV score", y="Density") +
  ggtitle("Public support") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  theme(legend.justification=c(1,1), legend.position=c(1,1)) + # put legend in top right corner
  theme(legend.text=element_text(size=8), legend.title=element_text(size=8))
lcv.yrange <- layer_scales(s1)$y$range$range[2]  # get y range of density plot
s1.errorbar <- s1 +   # add error bars to plot
  geom_errorbarh(aes(y=-lcv.yrange/10, xmax=lcv.mean.ppa + lcv.sd.ppa, xmin=lcv.mean.ppa - lcv.sd.ppa, height=lcv.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-lcv.yrange/10, x=lcv.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-lcv.yrange/13, xmax=lcv.mean.cpa + lcv.sd.ppa, xmin=lcv.mean.cpa - lcv.sd.cpa, height=lcv.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-lcv.yrange/13, x=lcv.mean.cpa), color="red", size=pointsize)

# Forestry sector
forest.vals <- get(paste0("mean.forestry.",bufdist1,"kmBuffer"))*100  # get mean and sd of value for each mode
forest.mean.ppa <- mean(forest.vals[which(DesMode=="President")], na.rm=TRUE)
forest.mean.cpa <- mean(forest.vals[which(DesMode=="Congress")], na.rm=TRUE)
forest.sd.ppa <- sd(forest.vals[which(DesMode=="President")], na.rm=TRUE)
forest.sd.cpa <- sd(forest.vals[which(DesMode=="Congress")], na.rm=TRUE)
s2 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0("mean.forestry.",bufdist1,"kmBuffer"))*100, fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x="Mean percentage of workforce", y="Density") +
  ggtitle("Forestry sector reliance") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
forest.yrange <- layer_scales(s2)$y$range$range[2]  # get y range of density plot
s2.errorbar <- s2 +   # add error bars to plot
  geom_errorbarh(aes(y=-forest.yrange/10, xmax=forest.mean.ppa + forest.sd.ppa, xmin=forest.mean.ppa - forest.sd.ppa, height=forest.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-forest.yrange/10, x=forest.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-forest.yrange/13, xmax=forest.mean.cpa + forest.sd.ppa, xmin=forest.mean.cpa - forest.sd.cpa, height=forest.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-forest.yrange/13, x=forest.mean.cpa), color="red", size=pointsize)

# Mining sector
mine.vals <- get(paste0("mean.mining.",bufdist1,"kmBuffer"))*100  # get mean and sd of value for each mode
mine.mean.ppa <- mean(mine.vals[which(DesMode=="President")], na.rm=TRUE)
mine.mean.cpa <- mean(mine.vals[which(DesMode=="Congress")], na.rm=TRUE)
mine.sd.ppa <- sd(mine.vals[which(DesMode=="President")], na.rm=TRUE)
mine.sd.cpa <- sd(mine.vals[which(DesMode=="Congress")], na.rm=TRUE)
s3 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0("mean.mining.",bufdist1,"kmBuffer"))*100, fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x="Mean percentage of workforce", y="Density") +
  ggtitle("Minerals sector reliance") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
mine.yrange <- layer_scales(s3)$y$range$range[2]  # get y range of density plot
s3.errorbar <- s3 +   # add error bars to plot
  geom_errorbarh(aes(y=-mine.yrange/10, xmax=mine.mean.ppa + mine.sd.ppa, xmin=mine.mean.ppa - mine.sd.ppa, height=mine.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-mine.yrange/10, x=mine.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-mine.yrange/13, xmax=mine.mean.cpa + mine.sd.ppa, xmin=mine.mean.cpa - mine.sd.cpa, height=mine.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-mine.yrange/13, x=mine.mean.cpa), color="red", size=pointsize)

#multiplot(s1, s2, s3, s4, cols=4)
multiplot(s1.errorbar, s2.errorbar, s3.errorbar, cols=1)



# MINIMUMS
bandwidth.multiplier <- 1   # adjust bandwidth (how much smoothing of histogram occurs - 1 is default)
pointsize <- 2   # adust size of points in error bars
heightfactor  <- 50  # adjust height of tick marks on end of error bars  (height = yrange of plot / heightfactor )
attach(PA_zonal.df)
bufdist1 <- "10"   # select buffer distance; must be character format, one of these choices: "10","20","50","100","250")
# LCV score
lcv.vals <- get(paste0("min.LCV.",bufdist1,"kmBuffer"))  # get mean and sd of value for each mode
lcv.mean.ppa <- mean(lcv.vals[which(DesMode=="President")], na.rm=TRUE)
lcv.mean.cpa <- mean(lcv.vals[which(DesMode=="Congress")], na.rm=TRUE)
lcv.sd.ppa <- sd(lcv.vals[which(DesMode=="President")], na.rm=TRUE)
lcv.sd.cpa <- sd(lcv.vals[which(DesMode=="Congress")], na.rm=TRUE)
s1 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0("min.LCV.",bufdist1,"kmBuffer")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x="Minimum LCV score", y="Density") +
  ggtitle("Public support") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  theme(legend.justification=c(1,1), legend.position=c(1,1)) + # put legend in top right corner
  theme(legend.text=element_text(size=8), legend.title=element_text(size=8))
lcv.yrange <- layer_scales(s1)$y$range$range[2]  # get y range of density plot
s1.errorbar <- s1 +   # add error bars to plot
  geom_errorbarh(aes(y=-lcv.yrange/10, xmax=lcv.mean.ppa + lcv.sd.ppa, xmin=lcv.mean.ppa - lcv.sd.ppa, height=lcv.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-lcv.yrange/10, x=lcv.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-lcv.yrange/13, xmax=lcv.mean.cpa + lcv.sd.ppa, xmin=lcv.mean.cpa - lcv.sd.cpa, height=lcv.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-lcv.yrange/13, x=lcv.mean.cpa), color="red", size=pointsize)

# Forestry sector
forest.vals <- get(paste0("min.forestry.",bufdist1,"kmBuffer"))*100  # get mean and sd of value for each mode
forest.mean.ppa <- mean(forest.vals[which(DesMode=="President")], na.rm=TRUE)
forest.mean.cpa <- mean(forest.vals[which(DesMode=="Congress")], na.rm=TRUE)
forest.sd.ppa <- sd(forest.vals[which(DesMode=="President")], na.rm=TRUE)
forest.sd.cpa <- sd(forest.vals[which(DesMode=="Congress")], na.rm=TRUE)
s2 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0("min.forestry.",bufdist1,"kmBuffer"))*100, fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x="Minimum percentage of workforce", y="Density") +
  ggtitle("Forestry sector reliance") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
forest.yrange <- layer_scales(s2)$y$range$range[2]  # get y range of density plot
s2.errorbar <- s2 +   # add error bars to plot
  geom_errorbarh(aes(y=-forest.yrange/10, xmax=forest.mean.ppa + forest.sd.ppa, xmin=forest.mean.ppa - forest.sd.ppa, height=forest.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-forest.yrange/10, x=forest.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-forest.yrange/13, xmax=forest.mean.cpa + forest.sd.ppa, xmin=forest.mean.cpa - forest.sd.cpa, height=forest.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-forest.yrange/13, x=forest.mean.cpa), color="red", size=pointsize)

# Mining sector
mine.vals <- get(paste0("min.mining.",bufdist1,"kmBuffer"))*100  # get mean and sd of value for each mode
mine.mean.ppa <- mean(mine.vals[which(DesMode=="President")], na.rm=TRUE)
mine.mean.cpa <- mean(mine.vals[which(DesMode=="Congress")], na.rm=TRUE)
mine.sd.ppa <- sd(mine.vals[which(DesMode=="President")], na.rm=TRUE)
mine.sd.cpa <- sd(mine.vals[which(DesMode=="Congress")], na.rm=TRUE)
s3 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0("min.mining.",bufdist1,"kmBuffer"))*100, fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x="Minimum percentage of workforce", y="Density") +
  ggtitle("Minerals sector reliance") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
mine.yrange <- layer_scales(s3)$y$range$range[2]  # get y range of density plot
s3.errorbar <- s3 +   # add error bars to plot
  geom_errorbarh(aes(y=-mine.yrange/10, xmax=mine.mean.ppa + mine.sd.ppa, xmin=mine.mean.ppa - mine.sd.ppa, height=mine.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-mine.yrange/10, x=mine.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-mine.yrange/13, xmax=mine.mean.cpa + mine.sd.ppa, xmin=mine.mean.cpa - mine.sd.cpa, height=mine.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-mine.yrange/13, x=mine.mean.cpa), color="red", size=pointsize)

#multiplot(s1, s2, s3, s4, cols=4)
multiplot(s1.errorbar, s2.errorbar, s3.errorbar, cols=1)



# MAXIMUMS
bandwidth.multiplier <- 1   # adjust bandwidth (how much smoothing of histogram occurs - 1 is default)
pointsize <- 2   # adust size of points in error bars
heightfactor  <- 50  # adjust height of tick marks on end of error bars  (height = yrange of plot / heightfactor )

attach(PA_zonal.df)
bufdist1 <- "10"   # select buffer distance; must be character format, one of these choices: "10","20","50","100","250")
# LCV score
lcv.vals <- get(paste0("max.LCV.",bufdist1,"kmBuffer"))  # get mean and sd of value for each mode
lcv.mean.ppa <- mean(lcv.vals[which(DesMode=="President")], na.rm=TRUE)
lcv.mean.cpa <- mean(lcv.vals[which(DesMode=="Congress")], na.rm=TRUE)
lcv.sd.ppa <- sd(lcv.vals[which(DesMode=="President")], na.rm=TRUE)
lcv.sd.cpa <- sd(lcv.vals[which(DesMode=="Congress")], na.rm=TRUE)
s1 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0("max.LCV.",bufdist1,"kmBuffer")) , fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x="Maximum LCV score", y="Density") +
  ggtitle("Public support") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  theme(legend.justification=c(1,1), legend.position=c(1,1)) + # put legend in top right corner
  theme(legend.text=element_text(size=8), legend.title=element_text(size=8))
lcv.yrange <- layer_scales(s1)$y$range$range[2]  # get y range of density plot
s1.errorbar <- s1 +   # add error bars to plot
  geom_errorbarh(aes(y=-lcv.yrange/10, xmax=lcv.mean.ppa + lcv.sd.ppa, xmin=lcv.mean.ppa - lcv.sd.ppa, height=lcv.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-lcv.yrange/10, x=lcv.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-lcv.yrange/13, xmax=lcv.mean.cpa + lcv.sd.ppa, xmin=lcv.mean.cpa - lcv.sd.cpa, height=lcv.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-lcv.yrange/13, x=lcv.mean.cpa), color="red", size=pointsize)

# Forestry sector
forest.vals <- get(paste0("max.forestry.",bufdist1,"kmBuffer"))*100  # get mean and sd of value for each mode
forest.mean.ppa <- mean(forest.vals[which(DesMode=="President")], na.rm=TRUE)
forest.mean.cpa <- mean(forest.vals[which(DesMode=="Congress")], na.rm=TRUE)
forest.sd.ppa <- sd(forest.vals[which(DesMode=="President")], na.rm=TRUE)
forest.sd.cpa <- sd(forest.vals[which(DesMode=="Congress")], na.rm=TRUE)
s2 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0("max.forestry.",bufdist1,"kmBuffer"))*100, fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x="Maximum proportion of workforce", y="Density") +
  ggtitle("Forestry sector reliance") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
forest.yrange <- layer_scales(s2)$y$range$range[2]  # get y range of density plot
s2.errorbar <- s2 +   # add error bars to plot
  geom_errorbarh(aes(y=-forest.yrange/10, xmax=forest.mean.ppa + forest.sd.ppa, xmin=forest.mean.ppa - forest.sd.ppa, height=forest.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-forest.yrange/10, x=forest.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-forest.yrange/13, xmax=forest.mean.cpa + forest.sd.ppa, xmin=forest.mean.cpa - forest.sd.cpa, height=forest.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-forest.yrange/13, x=forest.mean.cpa), color="red", size=pointsize)

# Mining sector
mine.vals <- get(paste0("max.mining.",bufdist1,"kmBuffer"))*100  # get mean and sd of value for each mode
mine.mean.ppa <- mean(mine.vals[which(DesMode=="President")], na.rm=TRUE)
mine.mean.cpa <- mean(mine.vals[which(DesMode=="Congress")], na.rm=TRUE)
mine.sd.ppa <- sd(mine.vals[which(DesMode=="President")], na.rm=TRUE)
mine.sd.cpa <- sd(mine.vals[which(DesMode=="Congress")], na.rm=TRUE)
s3 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0("max.mining.",bufdist1,"kmBuffer"))*100, fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x="Maximum proportion of workforce", y="Density") +
  ggtitle("Minerals sector reliance") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
mine.yrange <- layer_scales(s3)$y$range$range[2]  # get y range of density plot
s3.errorbar <- s3 +   # add error bars to plot
  geom_errorbarh(aes(y=-mine.yrange/10, xmax=mine.mean.ppa + mine.sd.ppa, xmin=mine.mean.ppa - mine.sd.ppa, height=mine.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-mine.yrange/10, x=mine.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-mine.yrange/13, xmax=mine.mean.cpa + mine.sd.ppa, xmin=mine.mean.cpa - mine.sd.cpa, height=mine.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-mine.yrange/13, x=mine.mean.cpa), color="red", size=pointsize)

#multiplot(s1, s2, s3, s4, cols=4)
multiplot(s1.errorbar, s2.errorbar, s3.errorbar, cols=1)





####################################################################################################################################################
#### DENSITY PLOTS OF ALL VARIABLES IN ONE FIGURE
####################################################################################################################################################

# ADD bandwidth adjustment here
bandwidth.multiplier <- 1   # adjust bandwidth (how much smoothing of histogram occurs - 1 is default)
pointsize <- 2   # adust size of points in error bars
heightfactor  <- 50  # adjust height of tick marks on end of error bars  (height = yrange of plot / heightfactor )
stat <- "mean"   # choose within-PA summary stat; "mean", "min", or "max"
statname <- "Mean"  # write out statistic name as you want it to appear in plot legend (full name, capitalized first letter:  "Mean", "Minimum", or "Maximum)
attach(PA_zonal.df)

# bird richness
bird.vals <- get(paste0(stat,".rich.bird"))  # get mean and sd of value for each mode
bird.mean.ppa <- mean(bird.vals[which(DesMode=="President")], na.rm=TRUE)
bird.mean.cpa <- mean(bird.vals[which(DesMode=="Congress")], na.rm=TRUE)
bird.sd.ppa <- sd(bird.vals[which(DesMode=="President")], na.rm=TRUE)
bird.sd.cpa <- sd(bird.vals[which(DesMode=="Congress")], na.rm=TRUE)
r1 <- ggplot() +   # make density plot
  geom_density(data=PA_zonal.df, aes(x=get(paste0(stat,".rich.bird")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x=paste0(statname, " # species"), y="Density") +
  ggtitle("(A) Bird richness") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
bird.yrange <- layer_scales(r1)$y$range$range[2]  # get y range of density plot
r1.errorbar <- r1 +   # add error bars to plot
  geom_errorbarh(aes(y=-bird.yrange/10, xmax=bird.mean.ppa + bird.sd.ppa, xmin=bird.mean.ppa - bird.sd.ppa, height=bird.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-bird.yrange/10, x=bird.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-bird.yrange/13, xmax=bird.mean.cpa + bird.sd.ppa, xmin=bird.mean.cpa - bird.sd.cpa, height=bird.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-bird.yrange/13, x=bird.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r1)$y$range$range[2]/20,layer_scales(r1)$y$range$range[2]/20))

# mammal richness
mammal.vals <- get(paste0(stat,".rich.mammal"))  # get mean and sd of value for each mode
mammal.mean.ppa <- mean(mammal.vals[which(DesMode=="President")], na.rm=TRUE)
mammal.mean.cpa <- mean(mammal.vals[which(DesMode=="Congress")], na.rm=TRUE)
mammal.sd.ppa <- sd(mammal.vals[which(DesMode=="President")], na.rm=TRUE)
mammal.sd.cpa <- sd(mammal.vals[which(DesMode=="Congress")], na.rm=TRUE)
r2 <- ggplot() +    
  geom_density(data=PA_zonal.df, aes(x=get(paste0(stat,".rich.mammal")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x=paste0(statname, " # species"), y="Density") +
  ggtitle("(B) Mammal richness") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
mammal.yrange <- layer_scales(r2)$y$range$range[2]  # get y range of density plot
r2.errorbar <- r2 +   # add error bars to plot
  geom_errorbarh(aes(y=-mammal.yrange/10, xmax=mammal.mean.ppa + mammal.sd.ppa, xmin=mammal.mean.ppa - mammal.sd.ppa, height=mammal.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-mammal.yrange/10, x=mammal.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-mammal.yrange/13, xmax=mammal.mean.cpa + mammal.sd.ppa, xmin=mammal.mean.cpa - mammal.sd.cpa, height=mammal.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-mammal.yrange/13, x=mammal.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r2)$y$range$range[2]/20,layer_scales(r2)$y$range$range[2]/20))

# fish richness
fish.vals <- get(paste0(stat,".rich.fish"))  # get mean and sd of value for each mode
fish.mean.ppa <- mean(fish.vals[which(DesMode=="President")], na.rm=TRUE)
fish.mean.cpa <- mean(fish.vals[which(DesMode=="Congress")], na.rm=TRUE)
fish.sd.ppa <- sd(fish.vals[which(DesMode=="President")], na.rm=TRUE)
fish.sd.cpa <- sd(fish.vals[which(DesMode=="Congress")], na.rm=TRUE)
r3 <- ggplot() +    
  geom_density(data=PA_zonal.df, aes(x=get(paste0(stat,".rich.fish")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x=paste0(statname, " # species"), y="Density") +
  ggtitle("(C) Fish richness") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  theme(legend.justification=c(1,1), legend.position=c(1,1)) + # put legend in top right corner
  theme(legend.text=element_text(size=8), legend.title=element_text(size=8))
fish.yrange <- layer_scales(r3)$y$range$range[2]  # get y range of density plot
r3.errorbar <- r3 +   # add error bars to plot
  geom_errorbarh(aes(y=-fish.yrange/10, xmax=fish.mean.ppa + fish.sd.ppa, xmin=fish.mean.ppa - fish.sd.ppa, height=fish.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-fish.yrange/10, x=fish.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-fish.yrange/13, xmax=fish.mean.cpa + fish.sd.ppa, xmin=fish.mean.cpa - fish.sd.cpa, height=fish.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-fish.yrange/13, x=fish.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r3)$y$range$range[2]/20,layer_scales(r3)$y$range$range[2]/20))

# amphibian richness
amphib.vals <- get(paste0(stat,".rich.amphib"))  # get mean and sd of value for each mode
amphib.mean.ppa <- mean(amphib.vals[which(DesMode=="President")], na.rm=TRUE)
amphib.mean.cpa <- mean(amphib.vals[which(DesMode=="Congress")], na.rm=TRUE)
amphib.sd.ppa <- sd(amphib.vals[which(DesMode=="President")], na.rm=TRUE)
amphib.sd.cpa <- sd(amphib.vals[which(DesMode=="Congress")], na.rm=TRUE)
r4 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0(stat,".rich.amphib")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x=paste0(statname, " # species"), y="Density") +
  ggtitle("(D) Amphibian richness") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
amphib.yrange <- layer_scales(r4)$y$range$range[2]  # get y range of density plot
r4.errorbar <- r4 +   # add error bars to plot
  geom_errorbarh(aes(y=-amphib.yrange/10, xmax=amphib.mean.ppa + amphib.sd.ppa, xmin=amphib.mean.ppa - amphib.sd.ppa, height=amphib.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-amphib.yrange/10, x=amphib.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-amphib.yrange/13, xmax=amphib.mean.cpa + amphib.sd.ppa, xmin=amphib.mean.cpa - amphib.sd.cpa, height=amphib.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-amphib.yrange/13, x=amphib.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r4)$y$range$range[2]/20,layer_scales(r4)$y$range$range[2]/20))

# reptile richness
reptile.vals <- get(paste0(stat,".rich.reptile"))  # get mean and sd of value for each mode
reptile.mean.ppa <- mean(reptile.vals[which(DesMode=="President")], na.rm=TRUE)
reptile.mean.cpa <- mean(reptile.vals[which(DesMode=="Congress")], na.rm=TRUE)
reptile.sd.ppa <- sd(reptile.vals[which(DesMode=="President")], na.rm=TRUE)
reptile.sd.cpa <- sd(reptile.vals[which(DesMode=="Congress")], na.rm=TRUE)
r5 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0(stat,".rich.reptile")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x=paste0(statname, " # species"), y="Density") +
  ggtitle("(E) Reptile richness") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
reptile.yrange <- layer_scales(r5)$y$range$range[2]  # get y range of density plot
r5.errorbar <- r5 +   # add error bars to plot
  geom_errorbarh(aes(y=-reptile.yrange/10, xmax=reptile.mean.ppa + reptile.sd.ppa, xmin=reptile.mean.ppa - reptile.sd.ppa, height=reptile.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-reptile.yrange/10, x=reptile.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-reptile.yrange/13, xmax=reptile.mean.cpa + reptile.sd.ppa, xmin=reptile.mean.cpa - reptile.sd.cpa, height=reptile.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-reptile.yrange/13, x=reptile.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r5)$y$range$range[2]/20,layer_scales(r5)$y$range$range[2]/20))

# tree richness
tree.vals <- get(paste0(stat,".rich.tree"))  # get mean and sd of value for each mode
tree.mean.ppa <- mean(tree.vals[which(DesMode=="President")], na.rm=TRUE)
tree.mean.cpa <- mean(tree.vals[which(DesMode=="Congress")], na.rm=TRUE)
tree.sd.ppa <- sd(tree.vals[which(DesMode=="President")], na.rm=TRUE)
tree.sd.cpa <- sd(tree.vals[which(DesMode=="Congress")], na.rm=TRUE)
r6 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0(stat,".rich.tree")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x=paste0(statname, " # species"), y="Density") +
  ggtitle("(F) Tree richness") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
tree.yrange <- layer_scales(r6)$y$range$range[2]  # get y range of density plot
r6.errorbar <- r6 +   # add error bars to plot
  geom_errorbarh(aes(y=-tree.yrange/10, xmax=tree.mean.ppa + tree.sd.ppa, xmin=tree.mean.ppa - tree.sd.ppa, height=tree.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-tree.yrange/10, x=tree.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-tree.yrange/13, xmax=tree.mean.cpa + tree.sd.ppa, xmin=tree.mean.cpa - tree.sd.cpa, height=tree.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-tree.yrange/13, x=tree.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r6)$y$range$range[2]/20,layer_scales(r6)$y$range$range[2]/20))

# G1/G2 richness
natserv.vals <- get(paste0(stat,".rich.natserv"))  # get mean and sd of value for each mode
natserv.mean.ppa <- mean(natserv.vals[which(DesMode=="President")], na.rm=TRUE)
natserv.mean.cpa <- mean(natserv.vals[which(DesMode=="Congress")], na.rm=TRUE)
natserv.sd.ppa <- sd(natserv.vals[which(DesMode=="President")], na.rm=TRUE)
natserv.sd.cpa <- sd(natserv.vals[which(DesMode=="Congress")], na.rm=TRUE)
r7 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0(stat,".rich.natserv")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x=paste0(statname, " rarity-weighted richness"), y="Density") +
  ggtitle("(G) G1 & G2 species richness") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
natserv.yrange <- layer_scales(r7)$y$range$range[2]  # get y range of density plot
r7.errorbar <- r7 +   # add error bars to plot
  geom_errorbarh(aes(y=-natserv.yrange/10, xmax=natserv.mean.ppa + natserv.sd.ppa, xmin=natserv.mean.ppa - natserv.sd.ppa, height=natserv.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-natserv.yrange/10, x=natserv.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-natserv.yrange/13, xmax=natserv.mean.cpa + natserv.sd.ppa, xmin=natserv.mean.cpa - natserv.sd.cpa, height=natserv.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-natserv.yrange/13, x=natserv.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r7)$y$range$range[2]/80,layer_scales(r7)$y$range$range[2]/80))

# ecological system richness
system.vals <- system.richness.rare  # get mean and sd of value for each mode
system.mean.ppa <- mean(system.vals[which(DesMode=="President")], na.rm=TRUE)
system.mean.cpa <- mean(system.vals[which(DesMode=="Congress")], na.rm=TRUE)
system.sd.ppa <- sd(system.vals[which(DesMode=="President")], na.rm=TRUE)
system.sd.cpa <- sd(system.vals[which(DesMode=="Congress")], na.rm=TRUE)
r8 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(system.richness.rare, fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x="Rarefied richness", y="Density") +
  ggtitle("(H) Ecological system richness") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
system.yrange <- layer_scales(r8)$y$range$range[2]  # get y range of density plot
r8.errorbar <- r8 +   # add error bars to plot
  geom_errorbarh(aes(y=-system.yrange/10, xmax=system.mean.ppa + system.sd.ppa, xmin=system.mean.ppa - system.sd.ppa, height=system.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-system.yrange/10, x=system.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-system.yrange/13, xmax=system.mean.cpa + system.sd.ppa, xmin=system.mean.cpa - system.sd.cpa, height=system.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-system.yrange/13, x=system.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r8)$y$range$range[2]/20,layer_scales(r8)$y$range$range[2]/20))

# climate refugial potential
climate.vals <- get(paste0(stat,".climate"))  # get mean and sd of value for each mode
climate.mean.ppa <- mean(climate.vals[which(DesMode=="President")], na.rm=TRUE)
climate.mean.cpa <- mean(climate.vals[which(DesMode=="Congress")], na.rm=TRUE)
climate.sd.ppa <- sd(climate.vals[which(DesMode=="President")], na.rm=TRUE)
climate.sd.cpa <- sd(climate.vals[which(DesMode=="Congress")], na.rm=TRUE)
r9 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0(stat,".climate")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x=paste0(statname, " refugial index"), y="Density") +
  ggtitle("(I) Climate refugial potential") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
climate.yrange <- layer_scales(r9)$y$range$range[2]  # get y range of density plot
r9.errorbar <- r9 +   # add error bars to plot
  geom_errorbarh(aes(y=-climate.yrange/10, xmax=climate.mean.ppa + climate.sd.ppa, xmin=climate.mean.ppa - climate.sd.ppa, height=climate.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-climate.yrange/10, x=climate.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-climate.yrange/13, xmax=climate.mean.cpa + climate.sd.ppa, xmin=climate.mean.cpa - climate.sd.cpa, height=climate.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-climate.yrange/13, x=climate.mean.cpa), color="red", size=pointsize) +
  scale_y_continuous(expand = c(layer_scales(r9)$y$range$range[2]/40,layer_scales(r9)$y$range$range[2]/40))

bufdist1 <- "10"   # select buffer distance; must be character format, one of these choices: "10","20","50","100","250")
# LCV score
lcv.vals <- get(paste0("mean.LCV.",bufdist1,"kmBuffer"))  # get mean and sd of value for each mode
lcv.mean.ppa <- mean(lcv.vals[which(DesMode=="President")], na.rm=TRUE)
lcv.mean.cpa <- mean(lcv.vals[which(DesMode=="Congress")], na.rm=TRUE)
lcv.sd.ppa <- sd(lcv.vals[which(DesMode=="President")], na.rm=TRUE)
lcv.sd.cpa <- sd(lcv.vals[which(DesMode=="Congress")], na.rm=TRUE)
s1 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0("mean.LCV.",bufdist1,"kmBuffer")), fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x="Mean LCV score", y="Density") +
  ggtitle("(J) Public support") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
lcv.yrange <- layer_scales(s1)$y$range$range[2]  # get y range of density plot
s1.errorbar <- s1 +   # add error bars to plot
  geom_errorbarh(aes(y=-lcv.yrange/10, xmax=lcv.mean.ppa + lcv.sd.ppa, xmin=lcv.mean.ppa - lcv.sd.ppa, height=lcv.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-lcv.yrange/10, x=lcv.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-lcv.yrange/13, xmax=lcv.mean.cpa + lcv.sd.ppa, xmin=lcv.mean.cpa - lcv.sd.cpa, height=lcv.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-lcv.yrange/13, x=lcv.mean.cpa), color="red", size=pointsize)

# Forestry sector
forest.vals <- get(paste0("mean.forestry.",bufdist1,"kmBuffer"))*100  # get mean and sd of value for each mode
forest.mean.ppa <- mean(forest.vals[which(DesMode=="President")], na.rm=TRUE)
forest.mean.cpa <- mean(forest.vals[which(DesMode=="Congress")], na.rm=TRUE)
forest.sd.ppa <- sd(forest.vals[which(DesMode=="President")], na.rm=TRUE)
forest.sd.cpa <- sd(forest.vals[which(DesMode=="Congress")], na.rm=TRUE)
s2 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0("mean.forestry.",bufdist1,"kmBuffer"))*100, fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x="Mean percentage of workforce", y="Density") +
  ggtitle("(K) Forestry sector reliance") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
forest.yrange <- layer_scales(s2)$y$range$range[2]  # get y range of density plot
s2.errorbar <- s2 +   # add error bars to plot
  geom_errorbarh(aes(y=-forest.yrange/10, xmax=forest.mean.ppa + forest.sd.ppa, xmin=forest.mean.ppa - forest.sd.ppa, height=forest.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-forest.yrange/10, x=forest.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-forest.yrange/13, xmax=forest.mean.cpa + forest.sd.ppa, xmin=forest.mean.cpa - forest.sd.cpa, height=forest.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-forest.yrange/13, x=forest.mean.cpa), color="red", size=pointsize)

# Mining sector
mine.vals <- get(paste0("mean.mining.",bufdist1,"kmBuffer"))*100  # get mean and sd of value for each mode
mine.mean.ppa <- mean(mine.vals[which(DesMode=="President")], na.rm=TRUE)
mine.mean.cpa <- mean(mine.vals[which(DesMode=="Congress")], na.rm=TRUE)
mine.sd.ppa <- sd(mine.vals[which(DesMode=="President")], na.rm=TRUE)
mine.sd.cpa <- sd(mine.vals[which(DesMode=="Congress")], na.rm=TRUE)
s3 <- ggplot() +
  geom_density(data=PA_zonal.df, aes(x=get(paste0("mean.mining.",bufdist1,"kmBuffer"))*100, fill=DesMode, color=DesMode), alpha=0.35, size=1, adjust=bandwidth.multiplier) + 
  labs(x="Mean percentage of workforce", y="Density") +
  ggtitle("(L) Minerals sector reliance") +
  scale_fill_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  scale_color_manual(name="Designation\nmode", breaks=c("Congress", "President"), values=c("red", "blue")) +
  theme(plot.title = element_text(size=11, face="bold")) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  guides(fill=FALSE, color=FALSE) # suppress legend
mine.yrange <- layer_scales(s3)$y$range$range[2]  # get y range of density plot
s3.errorbar <- s3 +   # add error bars to plot
  geom_errorbarh(aes(y=-mine.yrange/10, xmax=mine.mean.ppa + mine.sd.ppa, xmin=mine.mean.ppa - mine.sd.ppa, height=mine.yrange/heightfactor), color="blue", size=1) +
  geom_point(aes(y=-mine.yrange/10, x=mine.mean.ppa), color="blue", size=pointsize) +
  geom_errorbarh(aes(y=-mine.yrange/13, xmax=mine.mean.cpa + mine.sd.ppa, xmin=mine.mean.cpa - mine.sd.cpa, height=mine.yrange/heightfactor), color="red", size=1) +
  geom_point(aes(y=-mine.yrange/13, x=mine.mean.cpa), color="red", size=pointsize)


detach(PA_zonal.df)
multiplot(r1.errorbar, r2.errorbar, r3.errorbar,r4.errorbar,r5.errorbar,r6.errorbar,r7.errorbar,r8.errorbar,r9.errorbar, s1.errorbar, s2.errorbar, s3.errorbar, cols=3)   # density plots with means and SDs







################################################################################################
### MAPS OF ECOLOGICAL VARIABLES FOR SUPPLEMENTAL INFO
################################################################################################


# convert raster layers to dataframes to allow plotting in ggplot
rich.bird[rich.bird==0] <- NA
rich.bird.df <- as.data.frame(as(rich.bird, "SpatialPixelsDataFrame"))
colnames(rich.bird.df) <- c("value", "x", "y")
rich.mammal.df <- as.data.frame(as(rich.mammal, "SpatialPixelsDataFrame"))
colnames(rich.mammal.df) <- c("value", "x", "y")
rich.reptile.df <- as.data.frame(as(rich.reptile, "SpatialPixelsDataFrame"))
colnames(rich.reptile.df) <- c("value", "x", "y")
rich.tree.df <- as.data.frame(as(rich.tree, "SpatialPixelsDataFrame"))
colnames(rich.tree.df) <- c("value", "x", "y")
rich.natserv.df <- as.data.frame(as(rich.natserv, "SpatialPixelsDataFrame"))
colnames(rich.natserv.df) <- c("value", "x", "y")
climate.df <- as.data.frame(as(climate, "SpatialPixelsDataFrame"))
colnames(climate.df) <- c("value", "x", "y")

# create individual plots
birdplot <- ggplot() +
  geom_tile(data=rich.bird.df, aes(x=x, y=y, fill=value)) +
  theme(axis.text= element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        plot.title=element_text(hjust=0.5),
        panel.grid.major = element_blank()) +
  scale_fill_gradientn(colors = c("blue","plum1","red")) +
  ggtitle("Bird species richness") +
  theme(plot.title = element_text(size = 20)) +
  labs(fill="")

mammalplot <- ggplot() +
  geom_tile(data=rich.mammal.df, aes(x=x, y=y, fill=value)) +
  theme(axis.text= element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        plot.title=element_text(hjust=0.5),
        panel.grid.major = element_blank()) +
  scale_fill_gradientn(colors = c("blue", "plum1","red")) +
  ggtitle("Mammal species richness") +
  theme(plot.title = element_text(size = 20)) +
  labs(fill="")

fishplot <- ggplot() +
  geom_sf(data=rich.fish, aes(fill=Join_Count), color=NA) +
  theme(axis.text= element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        plot.title=element_text(hjust=0.5),
        panel.grid.major = element_line(color="white")) +
  scale_fill_gradientn(colors = c("blue","plum1","red")) +
  ggtitle("Fish species richness") +
  theme(plot.title = element_text(size = 20)) +
  labs(fill="")

amphibianplot <- ggplot() +
  geom_sf(data=rich.amphib, aes(fill=Join_Count), color=NA) +
  theme(axis.text= element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        plot.title=element_text(hjust=0.5),
        panel.grid.major = element_line(color="white")) +
  scale_fill_gradientn(colors = c("blue","plum1","red")) +
  ggtitle("Amphibian species richness") +
  theme(plot.title = element_text(size = 20)) +
  labs(fill="")

reptileplot <- ggplot() +
  geom_tile(data=rich.reptile.df, aes(x=x, y=y, fill=value)) +
  theme(axis.text= element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        plot.title=element_text(hjust=0.5),
        panel.grid.major = element_blank()) +
  scale_fill_gradientn(colors = c("blue","plum1","red")) +
  ggtitle("Reptile species richness") +
  theme(plot.title = element_text(size = 20)) +
  labs(fill="")

treeplot <- ggplot() +
  geom_tile(data=rich.tree.df, aes(x=x, y=y, fill=value)) +
  theme(axis.text= element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        plot.title=element_text(hjust=0.5),
        panel.grid.major = element_blank()) +
  scale_fill_gradientn(colors = c("blue","plum1","red")) +
  ggtitle("Tree species richness") +
  theme(plot.title = element_text(size = 20)) +
  labs(fill="")

G1G2plot <- ggplot() +
  geom_tile(data=rich.natserv.df, aes(x=x, y=y, fill=value)) +
  theme(axis.text= element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        plot.title=element_text(hjust=0.5),
        panel.grid.major = element_blank()) +
  scale_fill_gradientn(colors = c("blue","plum1","red")) +
  ggtitle("G1/G2 species richness") +
  theme(plot.title = element_text(size = 20)) +
  labs(fill="")

climateplot <- ggplot() +
  geom_tile(data=climate.df, aes(x=x, y=y, fill=value)) +
  theme(axis.text= element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        plot.title=element_text(hjust=0.5),
        panel.grid.major = element_blank()) +
  scale_fill_gradientn(colors = c("blue","plum1","red")) +
  ggtitle("Climate refugial potential") +
  theme(plot.title = element_text(size = 20)) +
  labs(fill="")  

grid.arrange(birdplot, mammalplot, fishplot, amphibianplot, reptileplot, treeplot, G1G2plot, climateplot, nrow=2)




################################################################################################
### MAPS OF SOCIOPOLITICAL VARIABLES FOR SUPPLEMENTAL INFO
################################################################################################


# forestry sector
forestry.stack <- stack(list.files("C:/Users/Tyler/Desktop/Monuments/GeneratedData/TimberRasters/", full.names=TRUE)[-1])  # read in annual forestry rasters as stacks
mean.forestry <- mean(forestry.stack, na.rm=TRUE)
pal <- colorRampPalette(c("blue","plum1","red"))
plot(mean.forestry*100, axes=FALSE, box=FALSE, col=pal(10), main="Forestry sector reliance")

# mining sector
mining.stack <- stack(list.files("C:/Users/Tyler/Desktop/Monuments/GeneratedData/MineRasters/MineRst/", full.names=TRUE)[-1])  # read in annual mining rasters as stack
mean.mining <- mean(mining.stack, na.rm=TRUE)
par(mar = c(2, 2, 2, 2))
plot(mean.mining*100, axes=FALSE, box=FALSE, col=pal(10), main="Minerals sector reliance")

# LCV score
lcv.names <- as.list(list.files("C:/Users/Tyler/Desktop/Monuments/GeneratedData/LCVRasters/", full.names=TRUE)[-c(1:16)])
lcv.rasters <- lapply(lcv.names, raster)
getMaxExtent <- function(rasters) { # function for getting maximum extent from list of rasters
  extents <- sapply(rasters, FUN = function(x) {
    raster::extent(x)
  })
  r <- raster(ext = extents[[1]], nrows = rasters[[1]]@nrows, ncols = rasters[[1]]@ncols)
  max_extent <- sapply(extents, FUN = function(x) {
    r <<- raster::extend(r, x)
  })
  raster::extent(r)
}
lcv.maxExtent <- getMaxExtent(lcv.rasters)  # get maximum extent
lcv.rasters.extend <- lapply(lcv.rasters, function(x) {extend(x, lcv.maxExtent, value=NA)}) # extend all rasters to maximum extent to allow stacking
LCV.stack <- stack(lcv.rasters.extend)  # stack mosaicked outputs with common extent
mean.LCV <- mean(LCV.stack, na.rm=TRUE)
plot(mean.LCV, axes=FALSE, box=FALSE, col=pal(10), main="Public support")




#########################################################################################
### PROPORTION MISSING DATA AS A FUNCTION OF BUFFER SIZE
#########################################################################################

# for sociopolitical variables
varnames <- c("LCV", "forestry", "mining")
propNA.mat <- matrix(NA, nrow=length(varnames), ncol=5)  # blank matrix to hold proportion NA for each variable and buffer distance
attach(PA_zonal.df)  
for(i in 1:length(varnames)){
  mean.names <- paste0("mean.",varnames[i],c(".10kmBuffer",".20kmBuffer",".50kmBuffer",".100kmBuffer",".250kmBuffer"))    # list of variable names for means
  for(j in 1:5){
    propNA.mat[i,j] <- sum(is.na(get(mean.names[j])))/length(get(mean.names[j]))
  }
}
detach(PA_zonal.df)  
propNA.df <- data.frame(bufferDist=c(10,20,50,100,250), LCVscore=propNA.mat[1,], forestry=propNA.mat[2,], mining=propNA.mat[3,], stringsAsFactors=FALSE)
propNA.gather <- gather(data=propNA.df, key=socVar, value=propNA, LCVscore:mining)
ggplot(data=propNA.gather, aes(x=bufferDist, y=propNA, color=socVar)) +
  geom_line(linetype=1) +
  geom_point() +
  scale_color_discrete(name="Variable",
                       breaks=c("LCVscore", "forestry", "mining"),
                       labels=c("LCV score", "% forestry", "% mineral extraction")) +
  labs(x="Buffer distance (km)", y="Proportion missing data")



# for ecological variables (no buffer, so a single value per variable)
ecovarnames <- c("mean.rich.bird","mean.rich.mammal", "mean.rich.fish", "mean.rich.amphib","mean.rich.reptile","mean.rich.tree","mean.rich.natserv","system.richness.rare","mean.climate")
attach(PA_zonal.df)
propNA.ecol <- rep(NA, length(ecovarnames))
for(i in 1:length(ecovarnames)){
  propNA.ecol[i] <- sum(is.na(get(ecovarnames[i])))/nrow(PA_zonal.df)
}
propNA.ecol.df <- data.frame(cbind(ecovarnames, round(propNA.ecol,3)), stringsAsFactors = FALSE)




#############################################################################################################################
### BOXPLOTS SHOWING SHOWING COMPARISONS OF SOCIOPOLITICAL VARIABLES AMONG DESIGNATION MODES FOR DIFFERENT BUFFER DISTANCES
#############################################################################################################################

varnames <- c("LCV","forestry","mining")  # choose variable (LCV, Farm, Forestry, Mine)
ylabnames <- c("LCV score", "Proportion workforce\nin forestry sector", "Proportion workforce\nin minerals extraction sector")

for(i in 1:length(varnames)){
  mean.varcols <- grep(paste0("mean.",varnames[i]), names(PA_zonal.df))  # column indices for mean values of selected variable
  mean.tidy <- PA_zonal.df %>%   # long-form version of dataframe for mean values
    gather(bufferKm, meanVal, mean.varcols, factor_key=TRUE)
  
  mean.buffer.plot <- ggplot() +  # boxplot of mean values as a function of buffer distance
    geom_boxplot(data=mean.tidy, aes(x=bufferKm, y=meanVal, fill=DesMode), width=0.5) +
    scale_x_discrete(labels=c("10","20","50","100","250")) +
    labs(y=ylabnames[i], x="Buffer distance (km)") +
    scale_color_discrete(name="Designation\nmode") +
    scale_fill_discrete(name="Designation\nmode") 
  assign(paste0("plot",i),mean.buffer.plot)
}

multiplot(plot1, plot2, plot3, cols=1)   # combine in single plot



#########################################################################################################################
### SUMMARY STATISTICS FOR EACH VARIABLE AND DESIGNATION CATEGORY
#########################################################################################################################

### Generate table of min, median, and max values of distributions for two designation categories
var.list <- c("clipped_area_ac","mean.rich.bird","mean.rich.mammal", "mean.rich.fish", "mean.rich.amphib","mean.rich.reptile","mean.rich.tree","mean.rich.natserv","system.richness.rare","mean.climate", "mean.LCV.10kmBuffer","mean.forestry.10kmBuffer","mean.mining.10kmBuffer")
min.mat <- matrix(NA, nrow=length(var.list), ncol=2)
mean.mat <- matrix(NA, nrow=length(var.list), ncol=2)
max.mat <- matrix(NA, nrow=length(var.list), ncol=2)
#desmode.list <- c("Congress","President","President then Congress")   
attach(PA_zonal.df)
for(i in 1:length(var.list)){
  min.mat[i,1] <- min(get(var.list[i])[which(DesMode=="Congress")], na.rm=TRUE)
  min.mat[i,2] <- min(get(var.list[i])[which(DesMode=="President")], na.rm=TRUE)
  mean.mat[i,1] <- mean(get(var.list[i])[which(DesMode=="Congress")], na.rm=TRUE)
  mean.mat[i,2] <- mean(get(var.list[i])[which(DesMode=="President")], na.rm=TRUE)
  max.mat[i,1] <- max(get(var.list[i])[which(DesMode=="Congress")], na.rm=TRUE)
  max.mat[i,2] <- max(get(var.list[i])[which(DesMode=="President")], na.rm=TRUE)
}
detach(PA_zonal.df)

output.df <- data.frame()
for(i in 1:nrow(min.mat)){
  output.df <- rbind(output.df,min.mat[i,],mean.mat[i,],max.mat[i,])
}
colnames(output.df) <- c("Congress","President")
row.names <- c()
for(j in 1:length(var.list)){
  row.names <- c(row.names, paste0(c("min.","mean.","max."),var.list[j]))
}
rownames(output.df) <- row.names
output.df.round <- round(output.df, 2)
write.table(output.df.round, "clipboard", sep="\t")



################################################################################################
### SUMMARY STATISTICS FOR MONUMENTS UNDER REVIEW
################################################################################################

# Convert PA raw values to percentiles (relative to all federal PAs) for each variable
var.list <- c("clipped_area_ac","mean.rich.bird","mean.rich.mammal", "mean.rich.fish", "mean.rich.amphib","mean.rich.reptile","mean.rich.tree","mean.rich.natserv","system.richness.rare","mean.climate", "mean.LCV.10kmBuffer","mean.forestry.10kmBuffer","mean.mining.10kmBuffer")
pct.mat <- matrix(NA, ncol=length(var.list), nrow=nrow(PA_zonal.df))
attach(PA_zonal.df)
for(i in 1:length(var.list)){
  newdata <- get(var.list[i])
  pct.mat[,i] <- round(rank(newdata, na.last="keep", ties.method="max")/length(newdata)*100, 1)
}
pct.df <- cbind(UnitName, InReview, as.data.frame(pct.mat))
names(pct.df) <- c("UnitName", "InReview", var.list)
detach(PA_zonal.df)
pct.df.review <- pct.df[which(pct.df$InReview=="Yes"),]  # percentiles for monuments under administrative review
write.table(pct.df.review, "clipboard", sep="\t")



