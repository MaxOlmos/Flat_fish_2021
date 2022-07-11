rm(list=ls())
# Libraries ---------------------------------------------------------------
library(tidyr)
library(tidyverse)
library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(viridis)
library(reshape2)
library(ggmap)
library("rnaturalearth")
library("rnaturalearthdata")
library(sf)
library(geosphere)
library(rgdal)



# Load shapefile ocean ----------------------------------------------------
# -------------------------------------------------------------------------

shape_path <- "C:/Users/Maxime/Documents/Git/Flat_fish_2021/01_data/ne_50m_ocean/"
coast_shapefile <- paste(shape_path, "ne_50m_ocean.shp", sep="")
ocean <- readOGR(coast_shapefile)
ocean_shapefile <- st_read(coast_shapefile)
ocean_sf <- st_transform(ocean_shapefile,crs=crs_ref)

#plot sf object
ggplot() + geom_sf(data=ocean_sf, colour = "black", fill = "blue")


shapefile <- st_read("02_transformed_data/EBSshelf/EBSshelf.shp")
EBS <- st_transform(shapefile,crs=crs_ref)

shapefile <- st_read("01_data/EBSslope/EBSslope.shp")
Slope <- st_transform(shapefile,crs=crs_ref)



# Define LL ---------------------------------------------------------------
# -------------------------------------------------------------------------

# - LL1 : global footprint ------------------------------------------------
file_output <- "02_transformed_data/observer"
load(paste0(file_output,"/CPUE_catchability_adfg.RData"))### Use this to draw points around your data

plot(CPUE_catchability$long, CPUE_catchability$lat,xlim=c(-180,-154),ylim=c(54,62),col="red")
plot(ocean,xlim=c(-180,-155),ylim=c(54,62),add=TRUE)
plot(Slope$geometry,add=T)
LL1 <- locator()
LL1$x[31]
LL1$y[32]<-55.1

# transform in sf
LL1_sf <- st_as_sf(as.data.frame(LL1),coords = c("x", "y")) 
st_crs(LL1_sf)<- crs_ref

LL1_polygon <- LL1_sf %>%  concaveman::concaveman(concavity = 2)

# plot polygon LL1
ggplot() + 
  geom_sf(data = LL1_polygon,col="blue") +
  geom_sf(data=world, col=NA, fill="black")+
  coord_sf(xlim = xlims, ylim = ylims,expand=F)#+
  #geom_point(aes(x=LL1$x[32],y=55.1),col="red")

# same crs
#st_crs(point)<- crs_ref
st_crs(LL1_polygon) <- crs_ref

# - LL2 : one point out of water ---------------------------------------------
file_output <- "02_transformed_data/observer"
load(paste0(file_output,"/CPUE_catchability_adfg.RData"))### Use this to draw points around your data

plot(CPUE_catchability$long, CPUE_catchability$lat,xlim=c(-180,-154),ylim=c(54,62),col="red")
plot(ocean,axes=F,xlim=c(-180,-155),ylim=c(54,62),add=TRUE)
plot(total_EBS$geometry,add=T)
LL2 <- locator()

# transform in sf
LL2_sf <- st_as_sf(as.data.frame(LL2),coords = c("x", "y")) 
st_crs(LL2_sf)<- crs_ref

LL2_polygon <- LL2_sf %>%  concaveman::concaveman(concavity = 2)

# plot polygon LL2
ggplot() + 
  geom_sf(data = LL2_polygon,col="blue") +
  geom_sf(data=world, col=NA, fill="black")+
  coord_sf(xlim = xlims, ylim = ylims,expand=F)




# -------------------------------------------------------------------------
# Intersection ------------------------------------------------------------
xlims <- c(-180,-156)#range(pretty(Data_Geostat$Lon))
ylims <- c(53,65)#range(pretty(Data_Geostat$Lat))

inters_LL1_bath <- st_difference(LL1_polygon,ocean_sf)
inters_LL2_bath <- st_union(inters_LL1_bath,LL2_polygon)

total_EBS_temp <- st_union(inters_LL2_bath,Slope)

ggplot() + 
  geom_sf(data = total_EBS_temp,col="blue",fill="blue",alpha=0.3)+
# geom_sf(data=inters_LL1_bath,fill="red",alpha=0.3)+
 #   geom_sf(data=Slope,fill=NA,alpha=0.1)+
  theme_bw(base_size = 8)


total_EBS <- st_union(EBS,total_EBS_temp)


ggplot() + 
  #geom_sf(data = inters_LL2_bath,col="blue",fill="blue")+
  geom_sf(data=world, col=NA, fill="black")+
   #geom_sf(data=EBS)+
  #geom_sf(data=Slope)+
  geom_sf(data=total_EBS,fill="lightblue",col=NA)+
  coord_sf(xlim = xlims, ylim = ylims,expand=F)+
  geom_point(data=CPUE_catchability,aes(x=long,y=lat),col="red")+
  theme_bw(base_size = 8)

## LL3 :Polygon cut
x_coord <- c(-179,  -179,-160,-160)
y_coord <- c(61.1,63,63,61.1)
xym <- cbind(x_coord, y_coord)
xym

p = Polygon(xym)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))
plot(sps)

LL3 <- st_as_sf(sps)
st_crs(LL3) <- crs_ref

total_EBS2 <-  st_difference(total_EBS,LL3)
#inters_LL3_bath <-  st_difference(inters_LL2_bath,LL3)

ggplot() + 
  #geom_sf(data=test)
 # geom_sf(data = inters_LL3_bath,col="blue",fill="blue")+
  # geom_sf(data=world, col=NA, fill="black")+
  #geom_sf(data=sps_sf)+
  # geom_sf(data=EBS)+
  #geom_sf(data=Slope)+
  geom_sf(data=total_EBS2,col="green")+
  coord_sf(xlim = xlims, ylim = ylims,expand=F)+
  geom_point(data=CPUE_catchability,aes(x=long,y=lat),col="red")+
  theme_bw(base_size = 8)


# - LL3 : one point out of water ---------------------------------------------
file_output <- "02_transformed_data/observer"
load(paste0(file_output,"/CPUE_catchability_adfg.RData"))### Use this to draw points around your data
plot(CPUE_catchability$long, CPUE_catchability$lat,xlim=c(-180,-154),ylim=c(54,62),col="red")
plot(ocean,axes=F,xlim=c(-180,-155),ylim=c(54,62),add=TRUE)
plot(total_EBS$geometry,add=T)
LL4 <- locator()

# transform in sf
LL4_sf <- st_as_sf(as.data.frame(LL4),coords = c("x", "y")) 
st_crs(LL4_sf)<- crs_ref

LL4_polygon <- LL4_sf %>%  concaveman::concaveman(concavity = 2)

grid <- st_union(total_EBS2,LL4_polygon)
save(grid, file="grid.Rdata")
ggplot() + 
  geom_sf(data=grid,col="green")+
  coord_sf(xlim = xlims, ylim = ylims,expand=F)+
  geom_point(data=CPUE_catchability,aes(x=long,y=lat),col="red")+
  theme_bw(base_size = 8)


grid_coord <- as.data.frame(st_coordinates(grid))
save(grid_coord, file="grid_coord.Rdata")

plot(grid_coord$X, grid_coord$Y)

##### L5
LL5 <- locator()

# transform in sf
LL5_sf <- st_as_sf(as.data.frame(LL5),coords = c("x", "y")) 
st_crs(LL5_sf)<- crs_ref
LL5_polygon <- LL5_sf %>%  concaveman::concaveman(concavity = 2)

grid2 <- st_union(grid,LL5_polygon)
ggplot() + 
  geom_sf(data=grid2,col="green")+
  coord_sf(xlim = xlims, ylim = ylims,expand=F)+
  geom_point(data=CPUE_catchability,aes(x=long,y=lat),col="red")+
  theme_bw(base_size = 8)


grid2_coord <- as.data.frame(st_coordinates(grid2))
plot(grid2_coord$X, grid2_coord$Y,type="l")


##### L6
LL6 <- locator()
# transform in sf
LL6_sf <- st_as_sf(as.data.frame(LL6),coords = c("x", "y")) 
st_crs(LL6_sf)<- crs_ref
LL6_polygon <- LL6_sf %>%  concaveman::concaveman(concavity = 2)

grid3 <- st_union(grid2,LL6_polygon)
ggplot() + 
  geom_sf(data=grid3,col="green")+
  coord_sf(xlim = xlims, ylim = ylims,expand=F)+
  geom_point(data=CPUE_catchability,aes(x=long,y=lat),col="red")+
  theme_bw(base_size = 8)


grid3_coord <- as.data.frame(st_coordinates(grid3))
plot(grid3_coord$X, grid3_coord$Y)

x <- grid3_coord %>% filter(Y<58,Y>57,X< -170)
plot(x$X,x$Y)


LL7 <- locator()
# transform in sf
LL7_sf <- st_as_sf(as.data.frame(LL7),coords = c("x", "y")) 
st_crs(LL7_sf)<- crs_ref
LL7_polygon <- LL7_sf %>%  concaveman::concaveman(concavity = 2)

grid4 <- st_union(grid3,LL7_polygon)
ggplot() + 
  geom_sf(data=grid4,col="green")+
  coord_sf(xlim = xlims, ylim = ylims,expand=F)+
  geom_point(data=CPUE_catchability,aes(x=long,y=lat),col="red")+
  theme_bw(base_size = 8)

grid4_coord <- as.data.frame(st_coordinates(grid4))
plot(grid4_coord$X, grid4_coord$Y)

x <- grid4_coord %>% filter(Y<58,Y>57,X< -170)
plot(x$X,x$Y)

save(grid4_coord, file="grid4_coord.Rdata")

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
region_extent <- data.frame(long=grid4_coord$X, lat=grid4_coord$Y)
str(region_extent)


### An example of how to create user-defined extrapolation
### regions (extents) for VAST.

### Cecilia O'Leary and Cole Monnahan | December 2020

### The extrapolation region defines the extent over which the
### model predictions are integrated. Any density outside of this
### region will not be included in estimates of the index. It is
### not used in model fitting. It comes with many built-in
### regions but often a user needs to define their own. Here, we
### demonstrate two ways of doing this: (1) From a set of points
### representing the outer extent of the region; (2) from an
### existing shape file.

## > 'data.frame':	42 obs. of  2 variables:
## $ long: num  -166 -166 -165 -165 -164 ...
## $ lat : num  53.9 54.1 54.2 54.6 55 ...

#### Turn it into a spatial polygon object
## Need to duplicate a point so that it is connected
region_extent <- rbind(region_extent, region_extent[1,])
## https://www.maths.lancs.ac.uk/~rowlings/Teaching/Sheffield2013/cheatsheet.html
poly <- Polygon(region_extent)
polys <- Polygons(list(poly), ID='all')
sps <- SpatialPolygons(list(polys))
## I think the F_AREA could be dropped here
sps <- SpatialPolygonsDataFrame(sps, data.frame(Id=factor('all'), F_AREA=1, row.names='all'))
proj4string(sps)<- CRS("+proj=longlat +datum=WGS84")
sps <- spTransform(sps, CRS("+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))
### Get UTM zone for conversion to UTM projection
## retrieves spatial bounding box from spatial data [,1] is
## longitude
lon <- sum(bbox(sps)[1,])/2
## convert decimal degrees to utm zone for average longitude, use
## for new CRS
utmzone <- floor((lon + 180)/6)+1
crs_LL <- CRS('+proj=longlat +ellps=WGS84 +no_defs')
sps@proj4string <- crs_LL


### --------------------------------------------------
### Create the VAST extroplation grid for method 1 and 2
## Convert the final in polygon to UTM
crs_UTM <- CRS(paste0("+proj=utm +zone=",utmzone," +ellps=WGS84 +datum=WGS84 +units=m +no_defs "))
region_polygon <- spTransform(sps, crs_UTM)

### Construct the extroplation grid for VAST using sf package
## Size of grid **in meters** (since working in UTM). Controls
## the resolution of the grid.
cell_size <- 2000
## This step is slow at high resolutions
region_grid <- st_make_grid(region_polygon, cellsize = cell_size, what = "centers")
## Convert region_grid to Spatial Points to SpatialPointsDataFrame
region_grid <- as(region_grid, "Spatial")
region_grid_sp <- as(region_grid, "SpatialPointsDataFrame")
## combine shapefile data (region_polygon) with Spatial Points
## (region_grid_spatial) & place in SpatialPointsDataFrame data
## (this provides you with your strata identifier (here called
## Id) in your data frame))
region_grid_sp@data <- over(region_grid, region_polygon)

## Convert back to lon/lat coordinates as that is what VAST uses
region_grid_LL <- as.data.frame(spTransform(region_grid_sp, crs_LL))
region_df <- with(region_grid_LL,
                  data.frame(Lon=coords.x1,
                             Lat=coords.x2, Id,
                             Area_km2=( (cell_size/1000)^2),
                             row=1:nrow(region_grid_LL)))
## Filter out the grid that does not overlap (outside extent)
region <- subset(region_df, !is.na(Id))
## This is the final file needed.
str(region)
## > 'data.frame':	106654 obs. of  5 variables:
##  $ Lon     : num  -166 -166 -166 -166 -166 ...
##  $ Lat     : num  53.9 53.9 54 53.9 53.9 ...
##  $ Id      : Factor w/ 1 level "all": 1 1 1 1 1 1 1 1 1 1 ...
##  $ Area_km2: num  4 4 4 4 4 4 4 4 4 4 ...
##  $ row     : int  401 402 975 976 977 978 1549 1550 1551 1552 ...

### Save it to be read in and passed to VAST later.
saveRDS(region, file = "user_region.rds")
### End of creating user extrapolation region object
### --------------------------------------------------

### Quick plots of the process for method 1
png('user_region.png', width=7, height=7, units='in', res=200)
par(mfrow=c(2,2))
with(region_extent, plot(long, lat, main='Extent in points in LL'))
plot(region_polygon, main='Polygon in UTM', axes=TRUE)
plot(region_grid, col=ifelse(is.na(region_df$Id), 'red', 'black'),
     axes=TRUE, main='Extrapolation area UTM')
with(region, plot(Lon, Lat, main='Extrapolation region in LL', pch='.'))
dev.off()


LL <- list()
LL$LL1 <- LL1
LL$LL2 <- LL2
LL$LL3 <- LL3
LL$LL4 <- LL4
LL$LL5 <- LL5
LL$LL6 <- LL6
LL$LL7 <- LL7

save(LL, file="LL_tot.Rdata")


user_region <- readRDS('user_region.rds')
str(user_region)
plot(user_region$Lon,user_region$Lon)
plot(user_region)
plot(user_region$Lon, user_region$Lat, main='Extrapolation region in LL', pch='.')
