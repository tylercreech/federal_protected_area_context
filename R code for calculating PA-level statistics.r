# National Monuments Analysis - Generate zonal statistics for protected areas


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
library(matrixStats)
library(velox)
library(doParallel)
library(foreach)
library(snow)

infolder <- "C:/Users/Tyler/Desktop/Monuments/GeneratedData"  # data folder


###################################################################################################################################################
# specify which geographic subset you want to work with:
subset <- "lower48"   # options are "lower48", "east", or "west"
###################################################################################################################################################



### READ IN PREPPED SPATIAL DATA
rich.bird <- raster(paste(infolder, "/rich.bird.tif", sep=""))  # bird richness
rich.mammal <- raster(paste(infolder, "/rich.mammal.tif", sep=""))  # mammal richness
rich.tree <- raster(paste(infolder, "/rich.tree.tif", sep=""))  # tree richness
rich.reptile <- raster(paste(infolder, "/rich.reptile.tif", sep=""))  # reptile richness
rich.fish <- st_read(paste(infolder, "/rich.fish.shp", sep=""), stringsAsFactors=FALSE)  # fish richness
rich.amphib <- st_read(paste(infolder, "/rich.amphib.shp", sep=""), stringsAsFactors=FALSE)  # amphibian richness
rich.natserv <- raster(paste(infolder, "/natserv.tif", sep=""))  # rarity weighted richness from NatureServ
natlandcover <- raster("C:/Users/Tyler/Desktop/Monuments/GeneratedData/gap_landfire_nationalterrestrialecosystems2011/gap_landfire_nationalterrestrialecosystems2011.tif")  # ecological systems
impervious <- raster(paste(infolder, "/impervious.tif", sep=""))  # percent impervious surface
climate <- raster(paste(infolder, "/climate.tif", sep=""))   # climate refugial potential
PA <- st_read(paste(infolder, "/post1996_PAs_", subset, ".shp", sep=""), stringsAsFactors=FALSE)   # protected area polygons
# rename attributes (st_write abbreviates these to match ESRI character limitations)
names(PA) <- c("UnitName", "InReview", "CurDesType","CurDesAuth", "OriDesAuth", "CurDesYear", "AntiqYear", "area_m2", "area_ac", "DesMode", "clipped_area_m2", "clipped_area_ac", "clipped_fraction", "geometry")
st_crs(PA) <- proj4string(natlandcover)
PA.sp <- as(PA, "Spatial") # convert sf polygon layer to a spatial layer first (required for extract function)



### GET MEAN PERCENT IMPERVIOUS SURFACE FOR EACH PA
mean.impervious <- as.vector(raster::extract(impervious, y=PA.sp, fun=mean, na.rm=TRUE, progress="window"))  # extract percent impervious pixels within each PA
cell.count <- raster::extract(natlandcover, y=PA.sp, na.rm=TRUE, progress="window")  # get number of natlandcover pixels associated with each PA (will tell me how many can be sampled in rarefaction for ecological systems richness)
cell.count.vect <- unlist(lapply(cell.count, length))

# remove PAs that have >1/3 impervious surface or are Redesignated PAs
impervious.discards <- c("First Ladies National Historic Site", "New Bedford Whaling National Historical Park", "Stonewall National Monument", # list of PAs with >
              "Carter G. Woodson National Historic Site", "Pullman National Monument", "Birmingham Civil Rights National Monument",
              "Belmont-Paul Women's Equality National Monument", "Little Rock Central High School National Historic Site",
              "President William Jefferson Clinton Birthplace Home National Historic Site", "Manhattan Project National Historical Park", 
              "African Burial Ground National Monument", "Ronald Reagan Boyhood Home National Historic Site")
desmode.discards <- c("Pinnacles National Park", "Great Sand Dunes National Park", "Great Sand Dunes National Preserve", "Black Canyon of the Gunnison National Park", 
                      "First State National Historical Park", "Harriet Tubman Underground Railroad National Historical Park")
discards <- c(impervious.discards, desmode.discards)
PA.keep <- PA[-which(PA$UnitName %in% discards),]
PA.discard <- PA[which(PA$UnitName %in% discards),]
PA.keep.sp <- as(PA.keep, "Spatial")


### PA ZONAL STATISTICS FOR RASTER INPUTS ###
inputnames <- c("climate", "rich.bird", "rich.mammal", "rich.tree", "rich.reptile", "rich.natserv")  # rasters for which we want to calculate zonal stats
for(i in 1:length(inputnames)) {  # calculate zonal stats for each input raster
  start <- Sys.time()
  zonalvals <- raster::extract(x=get(inputnames[i]), y=PA.keep.sp, weights=TRUE, progress="window")  # extract raster values and weights (e.g., cell area proportions) within each PA polygon
  prop.nonNA <- sapply(zonalvals, function(x) sum(x[which(is.na(x[,1])==FALSE),2]) / sum(x[,2]))  # get proportion of PA area that is non-NA (i.e., has value in raster layer)
  meanvals <- sapply(zonalvals, function(x) ifelse(sum(is.na(x[,1])==FALSE)==0, NA, weightedMean(x[,1], x[,2], na.rm=TRUE)))   # calculate weighted mean (excluding NA cells)
  meanvals[which(prop.nonNA<0.9)] <- NA   # set mean val to NA for PAs with data available for <90% of their area
  assign(paste("mean",inputnames[i],sep="."),meanvals)  # write to output variable
  maxvals <- sapply(zonalvals, function(x) ifelse(sum(is.na(x[,1])==FALSE)==0, NA, max(x[,1], na.rm=TRUE))) # calculate max
  maxvals[which(prop.nonNA<0.9)] <- NA    # set max val to NA for PAs with data available for <90% of their area
  assign(paste("max",inputnames[i],sep="."),maxvals)   # write to output variable
  minvals <- sapply(zonalvals, function(x) ifelse(sum(is.na(x[,1])==FALSE)==0, NA, min(x[,1], na.rm=TRUE)))# calculate min
  minvals[which(prop.nonNA<0.9)] <- NA    # set min val to NA for PAs with data available for <90% of their area
  assign(paste("min",inputnames[i],sep="."),minvals)   # write to output variable
  end <- Sys.time()
  process <- end - start   # calculate processing time
  print(paste0(i, "Process=", process))
}



### PA ZONAL STATISTICSS FOR ECOLOGICAL SYSTEM RICHNESS ###

# use rarefaction method to account for differences in PA area
ecol.systems <- raster::extract(natlandcover, PA.keep.sp, progress="window")  # list of ecological systems (by ID) within each PA
system.richness <- as.numeric(lapply(ecol.systems, function(x) length(unique(x[is.na(x)==FALSE])))) # calculate raw number of ecological systems (i.e., without controlling for PA area)
mincells <- min(unlist(lapply(ecol.systems, function(x) length(x[is.na(x)==FALSE]))))  # get smallest number of non-NA cells within a PA
nsamples <- 1000  # number of random samples you want to use for rarefaction
richness.mat <- matrix(nrow=length(ecol.systems), ncol=nsamples) # preallocate matrix to hold means of sample
for(i in 1:nsamples) {  # loop through 1000 random samples
  sample.data <- lapply(ecol.systems, function(x) sample(x[is.na(x)==FALSE], size=mincells, replace=FALSE))  # sample mincells (without replacement) from non-NA values for each PA
  sample.richness <- as.numeric(lapply(sample.data, function(x) length(unique(x)))) # calculate number of unique values in sample (i.e., richness)
  richness.mat[,i] <- sample.richness  # write richness values for sample i to matrix
}
system.richness.rare <- rowMeans(richness.mat)  # calculate mean across samples for each PA
# set richness value to NA for PAs with <90 percent non-NA data
prop.nonNA <- lapply(ecol.systems, function(x) 1 - sum(is.na(x)) / length(x))
system.richness[prop.nonNA<0.9] <- NA
system.richness.rare[prop.nonNA<0.9] <- NA



### PA ZONAL STATISTICS FOR VECTOR INPUTS ###

# Fish richness
temp.rich.fish <- st_intersection(rich.fish, PA.keep) %>%  # intersect PAs and richness polygons
  mutate(intersectPolyArea =  as.numeric(st_area(geometry))) %>%  # calculate areas of intersection polygons
  group_by(UnitName) %>%  # group intersection polygons by PA
  mutate(sumIntersectArea = sum(intersectPolyArea)) %>% # get the sum of intersect polygon areas associated with each PA (could be different than PA area because of missing data (e.g., watersheds with no fish richness))
  mutate(overlapProportion = intersectPolyArea/sumIntersectArea) %>% # get the proportion of the summed intersect areas associated with each intersect polygon (these are the "weights")
  group_by(UnitName) %>%
  summarise(weightedMean = weightedMean(Join_Count, overlapProportion, na.rm=TRUE), weightedMedian = weightedMedian(Join_Count, overlapProportion, na.rm=TRUE), max=max(Join_Count), min=min(Join_Count), prop.nonNA=mean(sumIntersectArea)/mean(clipped_area_m2))  # get weighted mean, weighted median, maximum, minimum, and proportion of the total PA area with non-NA values
# for those PAs with less than 90% coverage of non-NA richness data, assign overall NA value
temp.rich.fish$weightedMean[temp.rich.fish$prop.nonNA<0.9] <- NA
temp.rich.fish$weightedMedian[temp.rich.fish$prop.nonNA<0.9] <- NA
temp.rich.fish$max[temp.rich.fish$prop.nonNA<0.9] <- NA
temp.rich.fish$min[temp.rich.fish$prop.nonNA<0.9] <- NA
# deal with PAs that overlap blank spots in richness map, and are therefore not represented in results of above richness calculation
rich.fish.df <- data.frame(UnitName=temp.rich.fish$UnitName, weightedMean=temp.rich.fish$weightedMean, weightedMedian=temp.rich.fish$weightedMedian, max=temp.rich.fish$max, min=temp.rich.fish$min, stringsAsFactors=FALSE) # create dataframe out of temp.rich.fish
fish.PAnames <- rich.fish.df$UnitName  # get list of PA names in the richness output
all.PAnames <- PA.keep$UnitName  # get list of all PA names, including those missing from richness output
'%notin%' <- function(x,y) !(x %in% y)
missing.fish.PAnames <- PA.keep$UnitName[which(PA.keep$UnitName %notin% rich.fish.df$UnitName)]  # find PAs missing from richness output
rich.fish.df.corrected <- data.frame(UnitName=c(rich.fish.df$UnitName, missing.fish.PAnames), mean.rich.fish=c(rich.fish.df$weightedMean, rep(NA, length(missing.fish.PAnames))), median.rich.fish=c(rich.fish.df$weightedMedian, rep(NA, length(missing.fish.PAnames))), max.rich.fish=c(rich.fish.df$max, rep(NA, length(missing.fish.PAnames))), min.rich.fish=c(rich.fish.df$min, rep(NA, length(missing.fish.PAnames))))  # add missing PAs to new dataframe with NA for mean, min, and max richness value

# Amphibian richness
temp.rich.amphib <- st_intersection(rich.amphib, PA.keep) %>%  # intersect PAs and richness polygons
  mutate(intersectPolyArea =  as.numeric(st_area(geometry))) %>%  # calculate areas of intersection polygons
  group_by(UnitName) %>%  # group intersection polygons by PA
  mutate(sumIntersectArea = sum(intersectPolyArea)) %>% # get the sum of intersect polygon areas associated with each PA (could be different than PA area because of missing data (e.g., watersheds with no amphib richness))
  mutate(overlapProportion = intersectPolyArea/sumIntersectArea) %>% # get the proportion of the summed intersect areas associated with each intersect polygon (these are the "weights")
  group_by(UnitName) %>%
  summarise(weightedMean = weightedMean(Join_Count, overlapProportion, na.rm=TRUE), weightedMedian = weightedMedian(Join_Count, overlapProportion, na.rm=TRUE), max=max(Join_Count), min=min(Join_Count), prop.nonNA=mean(sumIntersectArea)/mean(clipped_area_m2))  # get weighted mean, weighted median, maximum, minimum, and proportion of the total PA area with non-NA values
# for those PAs with less than 90% coverage of non-NA richness data, assign overall NA value
temp.rich.amphib$weightedMean[temp.rich.amphib$prop.nonNA<0.9] <- NA
temp.rich.amphib$weightedMedian[temp.rich.amphib$prop.nonNA<0.9] <- NA
temp.rich.amphib$max[temp.rich.amphib$prop.nonNA<0.9] <- NA
temp.rich.amphib$min[temp.rich.amphib$prop.nonNA<0.9] <- NA
# deal with PAs that overlap blank spots in richness map, and are therefore not represented in results of above richness calculation
rich.amphib.df <- data.frame(UnitName=temp.rich.amphib$UnitName, weightedMean=temp.rich.amphib$weightedMean, weightedMedian=temp.rich.amphib$weightedMedian, max=temp.rich.amphib$max, min=temp.rich.amphib$min, stringsAsFactors=FALSE) # create dataframe out of temp.rich.amphib
amphib.PAnames <- rich.amphib.df$UnitName  # get list of PA names in the richness output
all.PAnames <- PA.keep$UnitName  # get list of all PA names, including those missing from richness output
'%notin%' <- function(x,y) !(x %in% y)
missing.amphib.PAnames <- PA.keep$UnitName[which(PA.keep$UnitName %notin% rich.amphib.df$UnitName)]  # find PAs missing from richness output
rich.amphib.df.corrected <- data.frame(UnitName=c(rich.amphib.df$UnitName, missing.amphib.PAnames), mean.rich.amphib=c(rich.amphib.df$weightedMean, rep(NA, length(missing.amphib.PAnames))), median.rich.amphib=c(rich.amphib.df$weightedMedian, rep(NA, length(missing.amphib.PAnames))), max.rich.amphib=c(rich.amphib.df$max, rep(NA, length(missing.amphib.PAnames))), min.rich.amphib=c(rich.amphib.df$min, rep(NA, length(missing.amphib.PAnames))))  # add missing PAs to new dataframe with NA for mean, min, and max richness value





# SOCIOPOLITICAL VARIABLE ZONAL STATS

# buffer PAs by difference distances (units are m)
PA.econproj <- st_transform(PA.keep, proj4string(forestry.stack)) # reproject PA layer to match economic layers
PA.buf10 <- st_buffer(PA.econproj, dist=10000)
PA.buf20 <- st_buffer(PA.econproj, dist=20000)
PA.buf50 <- st_buffer(PA.econproj, dist=50000)
PA.buf100 <- st_buffer(PA.econproj, dist=100000)
PA.buf250 <- st_buffer(PA.econproj, dist=250000)

# clip buffers to lower48
'%notin%' <- function(x,y) !(x %in% y)  #remove islands (HA, PR, VI) and AK, then dissolve
states <- st_read(paste(infolder,"/states2.shp",sep="")) %>%
  filter(STATE %notin% c("Hawaii","Puerto Rico","U.S. Virgin Islands", "Alaska"))
lower48 <- st_union(states, by_feature=FALSE)# create outline polygon for lower 48 states
lower48 <- st_transform(lower48, proj4string(forestry.stack))
croplayernames <- c("PA.buf10", "PA.buf20", "PA.buf50", "PA.buf100", "PA.buf250")  # names of layers you want to crop
lower48.sp <- as(lower48, "Spatial") # convert lower48 sf layer to sp (so extent can be extracted by crop function)
for(k in 1:length(croplayernames)) {
  sp.input <- as(get(croplayernames[k]), "Spatial")
  temp <- raster::crop(sp.input, lower48.sp, progress = 'text')
  assign(croplayernames[k], temp)
}
PA.buf.list <- list(PA.buf10, PA.buf20, PA.buf50, PA.buf100, PA.buf250) # combine buffer layers in a list
bufdists <- c("10km", "20km", "50km", "100km", "250km")

# Set up parallel processing
UseCores <- detectCores() - 1   # number of cores (leave one open for other processing tasks)
cl <- makeCluster(UseCores)  # initiate cluster
registerDoParallel(cl)   # register cluster


# FORESTRY: calculate mean, min, and max values for each PA (at each buffer level)
forestry.stack <- stack(list.files("C:/Users/Tyler/Desktop/Monuments/GeneratedData/TimberRasters/", full.names=TRUE))  # read in annual forestry rasters as stacks
forestry.out <- foreach(i=1:length(PA.buf.list), .packages=c("raster", "sf", "sp")) %dopar% {# loop through buffer distances (parallelize with foreach?)
  PA.buf.select <- PA.buf.list[[i]]  # pull the correct buffered PA layer
  out.df <- data.frame(mean.val=rep(NA,nrow(PA.buf.select)), min.val=rep(NA,nrow(PA.buf.select)), max.val=rep(NA,nrow(PA.buf.select)))
  for(j in 1:length(PA.buf.select)){ # loop through individual PAs
    PA.indiv.select <- PA.buf.select[j,] # pull the correct individual PA polygon
    min.year <- as.numeric(PA.indiv.select$CurDesYear) - 8  # start year for assessing econ for selected PA
    max.year <- min(2015,(as.numeric(PA.indiv.select$CurDesYear)-1))  # end year for assessing econ for select PA
      # I used min argument in above line, because we only have econ data through 2015, so code would otherwise break for PAs designated in 2017
    forestry.stack.subset <- subset(forestry.stack, subset=(min.year-1986):(max.year-1986))  # subset the stack to the appropriate set of 8 years
    bbox <- extent(as(PA.indiv.select, "Spatial"))  # get bounding box for polygon
    forestry.stack.subset.crop <- crop(forestry.stack.subset, bbox)  # crop stack to the minimum bounding box of the buffered polygon (will make following steps much faster)
    mean.rast <- calc(forestry.stack.subset.crop, fun=mean) # calculate pixel-wise mean across the years in the stack
    zonalvals <- unlist(raster::extract(x=mean.rast, y=PA.indiv.select))  # extract raster values and weights (e.g., cell area proportions) within each PA polygon
    prop.nonNA <- length(which(is.na(zonalvals)==FALSE)) / length(zonalvals)  # get proportion of PA area that is non-NA (i.e., has value in raster layer)
    if(prop.nonNA>=0.9) {  # for PAs with >90% data coverage, calculate statistics
      out.df$mean.val[j] <- ifelse(sum(is.na(zonalvals)==FALSE)==0, NA, mean(zonalvals, na.rm=TRUE)) 
      out.df$min.val[j] <- ifelse(sum(is.na(zonalvals)==FALSE)==0, NA, min(zonalvals, na.rm=TRUE))
      out.df$max.val[j] <- ifelse(sum(is.na(zonalvals)==FALSE)==0, NA, max(zonalvals, na.rm=TRUE))
    }
  }
  return(out.df)
}
# output is a list of dataframes (one per PA buffer layer), with each dataframe giving the mean, min, and max forestry values for cells within each PA
# convert dataframe to set of vectors with informative names
for(i in 1:length(forestry.out)){ # loop through dataframes (one per buffer distance)
  df <- forestry.out[[i]]
  assign(paste0("mean.forestry.",bufdists[i],"Buffer"), df$mean.val)
  assign(paste0("min.forestry.",bufdists[i],"Buffer"), df$min.val)
  assign(paste0("max.forestry.",bufdists[i],"Buffer"), df$max.val)
}


# MINING: calculate mean, min, and max values for each PA (at each buffer level)
mining.stack <- stack(list.files("C:/Users/Tyler/Desktop/Monuments/GeneratedData/MineRasters/MineRst/", full.names=TRUE))  # read in annual mining rasters as stack
mining.out <- foreach(i=1:length(PA.buf.list), .packages=c("raster", "sf", "sp")) %dopar% {# loop through buffer distances (parallelize with foreach?)
  PA.buf.select <- PA.buf.list[[i]]  # pull the correct buffered PA layer
  out.df <- data.frame(mean.val=rep(NA,nrow(PA.buf.select)), min.val=rep(NA,nrow(PA.buf.select)), max.val=rep(NA,nrow(PA.buf.select)))
  for(j in 1:length(PA.buf.select)){ # loop through individual PAs
    PA.indiv.select <- PA.buf.select[j,] # pull the correct individual PA polygon
    min.year <- as.numeric(PA.indiv.select$CurDesYear) - 8  # start year for assessing econ for selected PA
    max.year <- min(2015,(as.numeric(PA.indiv.select$CurDesYear)-1))  # end year for assessing econ for select PA
    # I used min argument in above line, because we only have econ data through 2015, so code would otherwise break for PAs designated in 2017
    mining.stack.subset <- subset(mining.stack, subset=(min.year-1986):(max.year-1986))  # subset the stack to the appropriate set of 8 years
    bbox <- extent(as(PA.indiv.select, "Spatial"))  # get bounding box for polygon
    mining.stack.subset.crop <- crop(mining.stack.subset, bbox)  # crop stack to the minimum bounding box of the buffered polygon (will make following steps much faster)
    mean.rast <- calc(mining.stack.subset.crop, fun=mean) # calculate pixel-wise mean across the years in the stack
    zonalvals <- unlist(raster::extract(x=mean.rast, y=PA.indiv.select))  # extract raster values and weights (e.g., cell area proportions) within each PA polygon
    prop.nonNA <- length(which(is.na(zonalvals)==FALSE)) / length(zonalvals)  # get proportion of PA area that is non-NA (i.e., has value in raster layer)
    if(prop.nonNA>=0.9) {  # for PAs with >90% data coverage, calculate statistics
      out.df$mean.val[j] <- ifelse(sum(is.na(zonalvals)==FALSE)==0, NA, mean(zonalvals, na.rm=TRUE))
      out.df$min.val[j] <- ifelse(sum(is.na(zonalvals)==FALSE)==0, NA, min(zonalvals, na.rm=TRUE))
      out.df$max.val[j] <- ifelse(sum(is.na(zonalvals)==FALSE)==0, NA, max(zonalvals, na.rm=TRUE))
    }
  }
  return(out.df)
}
# output is a list of dataframes (one per PA buffer layer), with each dataframe giving the mean, min, and max mining values for cells within each PA
# convert dataframe to set of vectors with informative names
for(i in 1:length(mining.out)){ # loop through dataframes (one per buffer distance)
  df <- mining.out[[i]]
  assign(paste0("mean.mining.",bufdists[i],"Buffer"), df$mean.val)
  assign(paste0("min.mining.",bufdists[i],"Buffer"), df$min.val)
  assign(paste0("max.mining.",bufdists[i],"Buffer"), df$max.val)
}


# LCV: calculate mean, min, and max values for each PA (at each buffer level)
# read in LCV rasters - need to change to have common extent before stacking
lcv.names <- as.list(list.files("C:/Users/Tyler/Desktop/Monuments/GeneratedData/LCVRasters/", full.names=TRUE))
lcv.rasters <- lapply(lcv.names, raster)
# write function for getting maximum extent from list of rasters (to be called within for loop)
getMaxExtent <- function(rasters) { 
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

LCV.out <- foreach(i=1:length(PA.buf.list), .packages=c("raster", "sf", "sp")) %do% {# parallel version maxes out HAL's memory, so use sequential version
#LCV.out <- foreach(i=1:length(PA.buf.list), .packages=c("raster", "sf", "sp")) %dopar% {# loop through buffer distances (parallelize with foreach?)
  PA.buf.select <- PA.buf.list[[i]]  # pull the correct buffered PA layer
  out.df <- data.frame(mean.val=rep(NA,nrow(PA.buf.select)), min.val=rep(NA,nrow(PA.buf.select)), max.val=rep(NA,nrow(PA.buf.select)))
  for(j in 1:length(PA.buf.select)){ # loop through individual PAs
    PA.indiv.select <- PA.buf.select[j,] # pull the correct individual PA polygon
    min.year <- as.numeric(PA.indiv.select$CurDesYear) - 8  # start year for assessing LCV for selected PA
    max.year <- min(2016,(as.numeric(PA.indiv.select$CurDesYear)-1))  # end year for assessing LCV for select PA
    # I used min argument in above line, because we only have LCV data through 2016, so code would otherwise break for PAs designated in 2017
    LCV.stack.subset <- subset(LCV.stack, subset=(min.year-1971):(max.year-1971))  # subset the stack to the appropriate set of 8 years
    bbox <- extent(as(PA.indiv.select, "Spatial"))  # get bounding box for polygon
    LCV.stack.subset.crop <- crop(LCV.stack.subset, bbox)  # crop stack to the minimum bounding box of the buffered polygon (will make following steps much faster)
    mean.rast <- calc(LCV.stack.subset.crop, fun=mean) # calculate pixel-wise mean across the years in the stack
    zonalvals <- unlist(raster::extract(x=mean.rast, y=PA.indiv.select))  # extract raster values and weights (e.g., cell area proportions) within each PA polygon
    prop.nonNA <- length(which(is.na(zonalvals)==FALSE)) / length(zonalvals)  # get proportion of PA area that is non-NA (i.e., has value in raster layer)
    if(prop.nonNA>=0.9) {  # for PAs with >90% data coverage, calculate statistics
      out.df$mean.val[j] <- ifelse(sum(is.na(zonalvals)==FALSE)==0, NA, mean(zonalvals, na.rm=TRUE)) 
      out.df$min.val[j] <- ifelse(sum(is.na(zonalvals)==FALSE)==0, NA, min(zonalvals, na.rm=TRUE))
      out.df$max.val[j] <- ifelse(sum(is.na(zonalvals)==FALSE)==0, NA, max(zonalvals, na.rm=TRUE))
    }
    rm(LCV.stack.subset, bbox, LCV.stack.subset.crop, zonalvals, PA.indiv.select)
    gc()
  }
  return(out.df)
}
# output is a list of dataframes (one per PA buffer layer), with each dataframe giving the mean, min, and max forestry values for cells within each PA
# convert dataframe to set of vectors with informative names
for(i in 1:length(LCV.out)){ # loop through dataframes (one per buffer distance)
  df <- LCV.out[[i]]
  assign(paste0("mean.LCV.",bufdists[i],"Buffer"), df$mean.val)
  assign(paste0("min.LCV.",bufdists[i],"Buffer"), df$min.val)
  assign(paste0("max.LCV.",bufdists[i],"Buffer"), df$max.val)
}


stopCluster(cl)


### REORDER OUTPUTS FROM RASTER OPERATIONS ALPHABETICALLY ###

# Output variables from raster extract are sorted by FID, but outputs from sp operations were sorted alphabetically because the dplyr function "group_by" was used and it automatically alphabetizes results
names.by.alpha <- sort(PA.keep$UnitName)  # alphabetical vector of unit names
names.by.fid <- PA.keep$UnitName   # vector of unit names by FID
reorder <- match(names.by.alpha, names.by.fid)  # order in which elements from raster extract outputs should appear to be alphabetically ordered
mean.climate <- mean.climate[reorder]
max.climate <- max.climate[reorder]
min.climate <- min.climate[reorder]
mean.rich.bird <- mean.rich.bird[reorder]
max.rich.bird <- max.rich.bird[reorder]
min.rich.bird <- min.rich.bird[reorder]
mean.rich.mammal <- mean.rich.mammal[reorder]
max.rich.mammal <- max.rich.mammal[reorder]
min.rich.mammal <- min.rich.mammal[reorder]
mean.rich.tree <- mean.rich.tree[reorder]
max.rich.tree <- max.rich.tree[reorder]
min.rich.tree <- min.rich.tree[reorder]
mean.rich.reptile <- mean.rich.reptile[reorder]
max.rich.reptile <- max.rich.reptile[reorder]
min.rich.reptile <- min.rich.reptile[reorder]
mean.rich.natserv <- mean.rich.natserv[reorder]
max.rich.natserv <- max.rich.natserv[reorder]
min.rich.natserv <- min.rich.natserv[reorder]
mean.impervious <- mean.impervious[reorder]
mean.forestry.10kmBuffer <- mean.forestry.10kmBuffer[reorder]
mean.forestry.20kmBuffer <- mean.forestry.20kmBuffer[reorder]
mean.forestry.50kmBuffer <- mean.forestry.50kmBuffer[reorder]
mean.forestry.100kmBuffer <- mean.forestry.100kmBuffer[reorder]
mean.forestry.250kmBuffer <- mean.forestry.250kmBuffer[reorder]
min.forestry.10kmBuffer <- min.forestry.10kmBuffer[reorder]
min.forestry.20kmBuffer <- min.forestry.20kmBuffer[reorder]
min.forestry.50kmBuffer <- min.forestry.50kmBuffer[reorder]
min.forestry.100kmBuffer <- min.forestry.100kmBuffer[reorder]
min.forestry.250kmBuffer <- min.forestry.250kmBuffer[reorder]
max.forestry.10kmBuffer <- max.forestry.10kmBuffer[reorder]
max.forestry.20kmBuffer <- max.forestry.20kmBuffer[reorder]
max.forestry.50kmBuffer <- max.forestry.50kmBuffer[reorder]
max.forestry.100kmBuffer <- max.forestry.100kmBuffer[reorder]
max.forestry.250kmBuffer <- max.forestry.250kmBuffer[reorder]
mean.mining.10kmBuffer <- mean.mining.10kmBuffer[reorder]
mean.mining.20kmBuffer <- mean.mining.20kmBuffer[reorder]
mean.mining.50kmBuffer <- mean.mining.50kmBuffer[reorder]
mean.mining.100kmBuffer <- mean.mining.100kmBuffer[reorder]
mean.mining.250kmBuffer <- mean.mining.250kmBuffer[reorder]
min.mining.10kmBuffer <- min.mining.10kmBuffer[reorder]
min.mining.20kmBuffer <- min.mining.20kmBuffer[reorder]
min.mining.50kmBuffer <- min.mining.50kmBuffer[reorder]
min.mining.100kmBuffer <- min.mining.100kmBuffer[reorder]
min.mining.250kmBuffer <- min.mining.250kmBuffer[reorder]
max.mining.10kmBuffer <- max.mining.10kmBuffer[reorder]
max.mining.20kmBuffer <- max.mining.20kmBuffer[reorder]
max.mining.50kmBuffer <- max.mining.50kmBuffer[reorder]
max.mining.100kmBuffer <- max.mining.100kmBuffer[reorder]
max.mining.250kmBuffer <- max.mining.250kmBuffer[reorder]
mean.LCV.10kmBuffer <- mean.LCV.10kmBuffer[reorder]
mean.LCV20kmBuffer <- mean.LCV.20kmBuffer[reorder]
mean.LCV.50kmBuffer <- mean.LCV.50kmBuffer[reorder]
mean.LCV.100kmBuffer <- mean.LCV.100kmBuffer[reorder]
mean.LCV.250kmBuffer <- mean.LCV.250kmBuffer[reorder]
min.LCV.10kmBuffer <- min.LCV.10kmBuffer[reorder]
min.LCV.20kmBuffer <- min.LCV.20kmBuffer[reorder]
min.LCV.50kmBuffer <- min.LCV.50kmBuffer[reorder]
min.LCV.100kmBuffer <- min.LCV.100kmBuffer[reorder]
min.LCV.250kmBuffer <- min.LCV.250kmBuffer[reorder]
max.LCV.10kmBuffer <- max.LCV.10kmBuffer[reorder]
max.LCV.20kmBuffer <- max.LCV.20kmBuffer[reorder]
max.LCV.50kmBuffer <- max.LCV.50kmBuffer[reorder]
max.LCV.100kmBuffer <- max.LCV.100kmBuffer[reorder]
max.LCV.250kmBuffer <- max.LCV.250kmBuffer[reorder]
system.richness <- system.richness[reorder]
system.richness.rare <- system.richness.rare[reorder]




### COMBINE OUTPUT VARIABLES IN A SINGLE DATAFRAME AND SAVE ###

PA.df <- tbl_df(PA.keep)[,-ncol(PA.keep)]  # convert to a tbl object (and strip out geometry field)
PA.df <- PA.df[order(PA.df$UnitName),]  # sort original dataframe alphabetically

outputvars <- c("mean.climate", "max.climate", "min.climate", "mean.rich.bird", "max.rich.bird", "min.rich.bird", "mean.rich.mammal", "max.rich.mammal", "min.rich.mammal", 
                "mean.rich.tree", "max.rich.tree", "min.rich.tree", "mean.rich.reptile", "max.rich.reptile", "min.rich.reptile", "mean.rich.natserv", "max.rich.natserv", 
                "min.rich.natserv", "system.richness", "system.richness.rare", "mean.impervious", "mean.forestry.10kmBuffer", "min.forestry.10kmBuffer", "max.forestry.10kmBuffer",
                "mean.forestry.20kmBuffer", "min.forestry.20kmBuffer", "max.forestry.20kmBuffer", "mean.forestry.50kmBuffer", "min.forestry.50kmBuffer", "max.forestry.50kmBuffer",
                "mean.forestry.100kmBuffer", "min.forestry.100kmBuffer", "max.forestry.100kmBuffer", "mean.forestry.250kmBuffer", "min.forestry.250kmBuffer", "max.forestry.250kmBuffer",
                "mean.mining.10kmBuffer", "min.mining.10kmBuffer", "max.mining.10kmBuffer", "mean.mining.20kmBuffer", "min.mining.20kmBuffer", "max.mining.20kmBuffer",
                "mean.mining.50kmBuffer", "min.mining.50kmBuffer", "max.mining.50kmBuffer","mean.mining.100kmBuffer", "min.mining.100kmBuffer", "max.mining.100kmBuffer",
                "mean.mining.250kmBuffer", "min.mining.250kmBuffer", "max.mining.250kmBuffer", "mean.LCV.10kmBuffer", "min.LCV.10kmBuffer", "max.LCV.10kmBuffer", "mean.LCV.20kmBuffer", "min.LCV.20kmBuffer", "max.LCV.20kmBuffer",
                "mean.LCV.50kmBuffer", "min.LCV.50kmBuffer", "max.LCV.50kmBuffer", "mean.LCV.100kmBuffer", "min.LCV.100kmBuffer", "max.LCV.100kmBuffer", "mean.LCV.250kmBuffer", "min.LCV.250kmBuffer", "max.LCV.250kmBuffer")  # vector of names of all output variables

for(i in 1:length(outputvars)){  # add each output variables as a new column in dataframe
  PA.df <- data.frame(PA.df, get(outputvars[i]))
}
names(PA.df)[(ncol(PA.df)-length(outputvars)+1):ncol(PA.df)] <- outputvars # give names to new output variables in dataframe
rich.fish.amphib.df <- merge(rich.fish.df.corrected, rich.amphib.df.corrected, by="UnitName")
PA.df <- merge(PA.df, rich.fish.amphib.df, by="UnitName")
PA_zonal.df <- PA.df
write.csv(PA_zonal.df, "C:/Users/Tyler/Desktop/PA_zonal_stats_post1996_lower48_28Oct2018.csv")  # write to csv
save(PA_zonal.df, file="C:/Users/Tyler/Desktop/PA_zonal_stats_post1996_lower48_28Oct2018.RData")  # output to workspace file
