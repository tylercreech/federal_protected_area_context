library(tidyverse)
library(tigris)
library(sf)
library(fasterize)

setwd("G:/My Drive/MonumentData/Generated Data/LCVScores/")
rt <- "http://scorecard.lcv.org/exports/"
yr <- seq(from=1972,to=2016,by=1)
hs <- "-house"
fl <- "-scorecard-grid-export.csv"

for(i in 1:length(yr)){
  y <- yr[i]
  link <- paste0(rt,y,hs,fl)
  fname <- paste0(y,hs,".csv")
  download.file(url=link, destfile=fname)
}

setwd("G:/My Drive/MonumentData/Generated Data/LCVScores/CongDistShp/")
rt <- "http://cdmaps.polisci.ucla.edu/shp/districts"
cn <- seq(from=92, to=114, by = 1)
cnp <- ifelse(cn < 100, paste0("0",cn), cn)
ext <- ".zip"

for(i in 1:length(cnp)){
  y <- cnp[i]
  link <- paste0(rt,y,ext)
  fname <- paste0("districts",y,ext)
  download.file(url=link, destfile=fname)
}

dist.zip <- list.files(getwd(), pattern=".zip")
lapply(dist.zip,unzip, junkpaths=TRUE)


##Load Cong District
setwd("G:/My Drive/MonumentData/Generated Data/LCVScores/CongDistShp/")

dist.files <- list.files(getwd(),pattern=".shp")
prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" #projection for NatLandCover in Monument Dataset

options(warn=-1)  # temporarily turn off warnings
for(i in 1:length(dist.files)) {
  temp1 <- readOGR(dist.files[i])
  temp1 <- spTransform(temp1, prj)
  if(length(which(gIsValid(temp1)==FALSE))==0) {
    print(paste("Geometry is valid for layer ",dist.files[i], sep=""))
    st_write(as(temp1,"sf"),paste0("G:/My Drive/MonumentData/Generated Data/LCVScores/CongDistShp/GeomClean/",dist.files[i]))
  } else {  # if invalid geometries are found (e.g., Ring Self-intersection), convert to sp and then add zero-width buffer
    print("Invalid geometry found - applying zero-width buffer...")
    temp2 <- gBuffer(temp1, byid=TRUE, width=0)  # add zero-width buffer
    if(length(which(gIsValid(temp2)==FALSE))==0) {  # check again for invalid geometries
      print(paste("Geometry corrected for layer ", dist.files[i], sep=""))
      temp3 <- as(temp2, "sf")
      st_write(temp3,paste0("G:/My Drive/MonumentData/Generated Data/LCVScores/CongDistShp/GeomClean/",dist.files[i]))
    } else {
      stop(paste("Unable to correct geometry for layer ",dist.files[i],"!!!", sep=""))
    }
    rm(temp1, temp2, temp3)
  }
}
options(warn=0)  # turn warnings back on

setwd("G:/My Drive/MonumentData/Generated Data/LCVScores/")
CongID <- sort(rep(cnp, 2))
CongID <- CongID[-1]
Cong.Lookup <- data.frame(CongID, yr, stringsAsFactors = FALSE)

LCVtext <- list.files("G:/My Drive/MonumentData/Generated Data/LCVScores/", pattern="house.csv")
excluded.states <- c("AS","GU","MP","PR","UM","VI","AK", "HI","DC")
exc.state <- subset(tigris::fips_codes, (fips_codes$state %in% excluded.states))
es <- unique(exc.state$state_name)

rast_LCV <- function(LCVFiles, CongLookup, excstate){
df1 <- read_csv(LCVFiles, skip=6, na = c("", "NA", "n/a")) %>% 
  mutate(st = substr(District, 1,2),
         dist = substr(District, 4,5),
         yr = substr(colnames(.)[4], 1,5)) %>% 
  dplyr::select(st, dist, yr, contains("Score")) %>% 
  dplyr::select(1:4)   
colnames(df1)[4] <- "LCVScore"  
df1$dist <- ifelse(df1$dist == "AL", "00", df1$dist)  

df <- df1 %>%
  group_by(st, dist,yr) %>% 
  summarise(LCVmn = mean(LCVScore, na.rm=TRUE)) %>% 
    left_join(., unique(fips_codes[,c(1,3)]), by = c("st"= "state")) 
  
con.lu <- CongLookup[CongLookup$yr == unique(as.numeric(df$yr)),1]

shp <- st_read(paste0("G:/My Drive/MonumentData/Generated Data/LCVScores/CongDistShp/GeomClean/districts",con.lu,".shp"), stringsAsFactors = FALSE) %>% 
  st_transform(., crs = prj)
shp$DISTRICT <- str_pad(shp$DISTRICT,2,"left","0")
LCV.shp <- shp %>% filter(.,!STATENAME %in% excstate)  %>% 
  left_join(., df, c("STATENAME" = "state_name", "DISTRICT"="dist"))
t.rst <- fasterize::raster(LCV.shp, res=270)
LCVrst <- fasterize::fasterize(LCV.shp, t.rst, field = "LCVmn", fun="min")
}

prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

LCVs1 <- map(LCVtext[1:10], function(x) rast_LCV(x, CongLookup = Cong.Lookup, excstate = es))
LCVsal <- map(LCVs1, function(x)raster::crop(raster::extend(LCVs1[[1]], x), LCVs1[[1]]))
LCVbrick1 <- do.call(raster::stack, LCVsal)
rm(LCVs1, LCVsal)

LCVs2 <- map(LCVtext[11:20], function(x) rast_LCV(x, CongLookup = Cong.Lookup, excstate = es))
LCVsal <- map(LCVs2, function(x)raster::crop(raster::extend(LCVs2[[1]], x), LCVs2[[1]]))
LCVbrick2 <- do.call(raster::stack, LCVsal)
rm(LCVs2, LCVsal)
LCV1.2 <- raster::stack(LCVbrick1, LCVbrick2)
rm(LCVbrick1,LCVbrick2)
names(LCV1.2) <- yr[1:20]
raster::writeRaster(LCV1.2, filename = paste0("G:/My Drive/MonumentData/Generated Data/LCVScores/Rasters/", names(LCV1.2)),bylayer=TRUE, format = "GTiff")
rm(LCV1.2)

gc()

LCVs3 <- map(LCVtext[21:30], function(x) rast_LCV(x, CongLookup = Cong.Lookup, excstate = es))
LCVsal <- map(LCVs3, function(x)raster::crop(raster::extend(LCVs3[[1]], x), LCVs3[[1]]))
LCVbrick3 <- do.call(raster::stack, LCVsal)
rm(LCVs3, LCVsal)

LCVs4 <- map(LCVtext[31:40], function(x) rast_LCV(x, CongLookup = Cong.Lookup, excstate = es))
LCVsal <- map(LCVs4, function(x)raster::crop(raster::extend(LCVs4[[1]], x), LCVs4[[1]]))
LCVbrick4 <- do.call(raster::stack, LCVsal)
LCV3al <- raster::crop(raster::extend(LCVbrick4,LCVbrick3), LCVbrick4)
rm(LCV3al, LCVbrick3)
LCV3.4 <- raster::stack(LCV3al, LCVbrick4)
names(LCV3.4) <- yr[21:40]
raster::writeRaster(LCV3.4, filename = paste0("G:/My Drive/MonumentData/Generated Data/LCVScores/Rasters/", names(LCV3.4)),bylayer=TRUE, format = "GTiff")
rm(LCV3.4)

LCVs5 <- map(LCVtext[41:45], function(x) rast_LCV(x, CongLookup = Cong.Lookup, excstate = es))
LCVsal <- map(LCVs5, function(x)raster::crop(raster::extend(LCVs5[[1]], x), LCVs5[[1]]))
LCV5st <- do.call(raster::stack, LCVsal)
names(LCV5st) <- yr[41:45]
raster::writeRaster(LCV5st, filename = paste0("G:/My Drive/MonumentData/Generated Data/LCVScores/Rasters/", names(LCV5st)),bylayer=TRUE, format = "GTiff")




