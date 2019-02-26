# National Monuments Analysis - Data Preparation for Ecological Variable Layers


library(raster)
library(sp)
library(sf)
library(rgeos)
library(rgdal)
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyverse)
library(rvest)
library(maptools)

infolder <- "C:/Users/Tyler/Google Drive/MonumentData/" # folder where spatial data input layers are stored;


######################################################################################################
### DATA PREP
######################################################################################################


#### LOAD DATA ###

# Federal lands (regardless of protection status, from PADUS-CBI)
fedlands <- readOGR(paste(infolder,"CBI_federal_lands.shp", sep=""))
fedlands <- as(fedlands, "sf")  # reading in directly as sf produces an error due to null geometries; read in with rgdal, then convert to sf
fedlands <- st_union(fedlands, by_feature=FALSE)
st_write(fedlands, "D:/Data/MonumentData/Generated Data/fedlands.shp")
fedlands <- st_read("D:/Data/MonumentData/Generated Data/fedlands.shp") # necessary to reload files

# backwards climate velocity (AdaptWest)
climate <- raster(paste(infolder,"bwvelocityrefugiaindex.asc", sep=""))  # climate refugial potential raster
crs(climate) <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"  # this layer is missing projection info, but should be in Lambert Azimuthal Equal Area

# species richness (Jenkins et al - note that data are vector format for fish and amphibians, raster for the rest)
rich.amphib <- st_read(paste(infolder,"Amphibians_total_richness.shp",sep=""))
rich.fish <- st_read(paste(infolder,"Fish_total_richness.shp",sep=""))
rich.bird <- raster(paste(infolder,"Birds_total_richness.tif",sep=""))
rich.mammal <- raster(paste(infolder,"Mammals_total_richness.tif",sep=""))
rich.reptile <- raster(paste(infolder,"Reptiles_total_richness.tif",sep=""))
rich.tree <- raster(paste(infolder,"Trees_total_richness.tif",sep=""))

# GAP land cover (ecological systems)
landcover <- raster(paste(infolder,"lwr48_v2_2_1", sep=""))
# remove open water, developed, ag, disturbed cover types
lc.att <- as.data.frame(levels(landcover))  # get attribute table (and convert to df)
IDs.to.remove <- lc.att$ID[which(lc.att$NVC_CLASS %in% c("Agricultural Vegetation","Developed & Other Human Use", "Open Water", "Recently Disturbed or Modified"))]
rcl <- matrix(nrow=length(IDs.to.remove), ncol=2)
rcl[,1] <- IDs.to.remove
rcl[,2] <- NA
rcl.df <- as.data.frame(rcl)
colnames(rcl.df) <- c("from", "to")
natlandcover <- subs(landcover, rcl.df, by=1, which=2, subsWithNA=FALSE)
writeRaster(natlandcover, "D:/Data/MonumentData/Generated Data/natlandcvr.tif",format="GTiff",prj=TRUE)
natlandcover <- raster("D:/Data/MonumentData/Generated Data/natlandcvr.tif") #necessary to reload files

# US states
'%notin%' <- function(x,y) !(x %in% y)  #remove islands (HA, PR, VI) and AK, then dissolve
states <- st_read(paste(infolder,"states_albers.shp",sep="")) %>%
  filter(STATE %notin% c("Hawaii","Puerto Rico","U.S. Virgin Islands", "Alaska"))
lower48 <- st_union(states, by_feature=FALSE)# create outline polygon for lower 48 states


### REPROJECT SPATIAL LAYERS TO COMMON PROJECTION ###

newproj <- proj4string(natlandcover)
shapefiles <- c("PA","lower48","rich.fish","rich.amphib","bailey","fedlands")  # names of shapefiles to reproject ,"sectordom","lcv"
cat.rasters <- c("natlandcover")  # names of categorical rasters to reproject
cont.rasters <- c("rich.mammal", "rich.bird", "rich.tree", "rich.reptile", "climate")   # names of continuous rasters to reproject

for(i in 1:length(shapefiles)) {   # reproject shapefiles
  reproj <- st_transform(get(shapefiles[i]), crs=newproj)
  assign(shapefiles[i], reproj)
}
for(j in 1:length(cat.rasters)) {   # reproject categorical rasters
  reproj2 <- projectRaster(get(cat.rasters[j]), crs=newproj, method="ngb")
  assign(cat.rasters[j], reproj2)
}
for(k in 1:length(cont.rasters)) {   # reproject categorical rasters
  reproj3 <- projectRaster(get(cont.rasters[k]), crs=newproj, method="bilinear")
  assign(cont.rasters[k], reproj3)
}


### CHECK FOR INVALID GEOMETRIES IN SF LAYERS ###

geomlayernames <- c("lower48","rich.fish","rich.amphib","bailey","fedlands")  # names of layers you want to check - Need to re-do Fedlands on its own
options(warn=-1)  # temporarily turn off warnings
for(i in 1:length(geomlayernames)) {
  if(length(which(st_is_valid(get(geomlayernames[i]))==FALSE))==0) {
    print(paste("Geometry is valid for layer ",geomlayernames[i], sep=""))
  } else {  # if invalid geometries are found (e.g., Ring Self-intersection), convert to sp and then add zero-width buffer
    print("Invalid geometry found - applying zero-width buffer...")
    temp1 <- as(get(geomlayernames[i]), "Spatial")  # convert to sp
    temp2 <- gBuffer(temp1, byid=TRUE, width=0)  # add zero-width buffer
    temp3 <- as(temp2, "sf")   # convert back to sf
    if(length(which(st_is_valid(temp3)==FALSE))==0) {  # check again for invalid geometries
      assign(geomlayernames[i], temp3)
      print(paste("Geometry corrected for layer ", geomlayernames[i], sep=""))
    } else {
      stop(paste("Unable to correct geometry for layer ",geomlayernames[i],"!!!", sep=""))
    }
    rm(temp1, temp2, temp3)
  }
}
options(warn=0)  # turn warnings back on


### CROP INPUT LAYERS TO LOWER 48 ###
rasterOptions(tmpdir = "D:/RastTemp", progress = "text", maxmemory = 2e+08, todisk=TRUE)
croplayernames <- c("rich.fish","rich.amphib","bailey" ,"fedlands")  # names of layers you want to crop
lower48.sp <- as(lower48, "Spatial") # convert lower48 sf layer to sp (so extent can be extracted by crop function)
for(k in 1:length(croplayernames)) {
  sp.input <- as(get(croplayernames[k]), "Spatial")
  temp4 <- raster::crop(sp.input, lower48.sp, progress = 'text')
  assign(croplayernames[k], temp4)
}

croprastnames <- c("natlandcover", "rich.mammal", "rich.bird", "rich.tree", "rich.reptile", "climate")
for(k in 1:length(croprastnames)) {
  sp.input <- get(croprastnames[k])
  temp4 <- raster::crop(sp.input, lower48.sp, progress = 'text')
  assign(croprastnames[k], temp4)
}


### WRITE CROPPED LAYERS TO FILE
for(l in 1:length(croplayernames)) {
  temp <- as(get(croplayernames[l]),"sf")
  st_write(temp, paste0("D:/Data/MonumentData/Generated Data/",croplayernames[l],".shp"))
}

for(l in 1:length(croprastnames)) {
  temp <- get(croprastnames[l])
  writeRaster(temp, paste0("D:/Data/MonumentData/Generated Data/",croprastnames[l],".tif"), format="GTiff", prj=TRUE)
}


### SAVE FINAL OBJECTS TO WORKSPACE
keepers <- c("lower48","rich.fish","rich.amphib","bailey","fedlands","sectordom","lcv","natlandcover","rich.mammal", "rich.bird", "rich.tree", "rich.reptile","climate")   # list of objects we want to save (this will be input to zonal stats script)
outfile <- "C:/Users/Tyler/Desktop/NatMon_prepped_data.RData"  # path for saved workspace
save(list=keepers, file=outfile)
