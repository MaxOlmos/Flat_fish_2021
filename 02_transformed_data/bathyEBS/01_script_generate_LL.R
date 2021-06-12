library(metR)
library(ggplot2)
library(VAST)                           # 3.8.0
library(ggplot2)                        # 2.10.0
library(dplyr)
library(tidyr)
theme_set(theme_bw())
library(concaveman)
library(sf)
library(sp) # 1.4.4
library(sf) # 0.9.6
getwd()

shape_path <- "C:/Users/Maxime/Documents/Git/Flat_fish/Pheno-flatfish/ADFGv2/shapefile/"
coast_shapefile <- paste(shape_path, "ne_50m_land.shp", sep="")
ocean <- readOGR(coast_shapefile)

bathy.dat<-read.table('C:/Users/Maxime/Documents/Git/Flat_fish_2021/01_data/bathyEBS/BeringDepth.txt',sep='') 
names(bathy.dat)<-c('lon','lat','depth')
bathy.dat$depth[bathy.dat$depth>0]<-NA#Avoid points above water
head(bathy.dat)
bathy.mat<- matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]


### Load data
file_output <- "02_transformed_data/observer"
load(paste0(file_output,"/CPUE_catchability_adfg.RData"))
### Use this to draw points around your data

plot(CPUE_catchability$long, CPUE_catchability$lat)
plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(54,62),add=TRUE)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(200,100,50),labcex=0.4,col='black',add=T)

LL <- locator()
LL2 <- cbind(LL$x, LL$y)
LL3 <- as.data.frame(LL2)
save(LL3,file="02_transformed_data/bathyEBS/LL3.RData")


polygon_bathy <- st_as_sf(as.data.frame(LL3),coords = c("V1", "V2"))%>% 
  concaveman::concaveman(concavity = 5/4)
point <- st_as_sf(bathy.dat,coords = c("lon","lat"))

bathy.dat2 <- st_intersection(point,polygon_bathy$polygons )  
bathy.dat3 <- cbind(st_coordinates(bathy.dat2),bathy.dat2$depth)
bathy.dat3 <- as.data.frame(bathy.dat3)
colnames(bathy.dat3) <- c("X","Y","depth")

ggplot()+ geom_contour( data=bathy.dat3,aes(x = X, y = Y,z = depth,group=NULL),breaks  = c(-200,-100,-50),col="black",size=1)
contour(unique(bathy.dat3$X),sort(unique(bathy.dat$Y)),bathy.mat,levels=-c(200,100,50),labcex=0.4,col='black',add=T)
