rm(list=ls())

library(TMB)
library(dplyr)
library(VAST)
library(FishStatsUtils)
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

# Shapefile ---------------------------------------------------------------
# -------------------------------------------------------------------------

# world
world <- ne_countries(scale = "medium", returnclass = "sf")
crs_ref <- "+proj=longlat +ellps=WGS84  +datum=WGS84 +no_defs"
crs_ref_utm <- "+proj=utm +ellps=WGS84  +units=km +zone=2 +datum=WGS84 +no_defs"
world <- st_transform(world,crs=crs_ref)
getwd()
# EBS
shapefile <- st_read("01_data/EBSshelf/EBSshelf.shp")
EBS <- st_transform(shapefile,crs=crs_ref)

xlims <- c(-180,-156)#range(pretty(Data_Geostat$Lon))
ylims <- c(54.5,65)#range(pretty(Data_Geostat$Lat))

# EBS
shapefile <- st_read("01_data/EBSslope/EBSslope.shp")
Slope <- st_transform(shapefile,crs=crs_ref)

xlims <- c(-180,-156)#range(pretty(Data_Geostat$Lon))
ylims <- c(53.5,65)#range(pretty(Data_Geostat$Lat))


p <- ggplot() + 

  geom_sf(data=world, col=NA, fill="black")+
  geom_sf(data=Slope,col="orange",fill=NA)+
  geom_sf(data = EBS,col="red",alpha=0.5,fill=NA)+
  coord_sf(xlim = xlims, ylim = ylims,expand=F)+
  theme_bw(base_size = 8)
p
