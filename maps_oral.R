# Load packages -----------------------------------------------------------
# -------------------------------------------------------------------------
library(VAST)
library(TMB)
library(dplyr)
library(rgdal)
library(raster)
library(pander)
library(splines)
#library(FishStatsUtils)
library(devtools)
library(sf)
theme_set(theme_bw())


# Maps --------------------------------------------------------------------
# -------------------------------------------------------------------------
library("rnaturalearth")
library("rnaturalearthdata")
library("scales")
library(ggplot2)

getwd()
world <- ne_countries(scale = "medium", returnclass = "sf")
shapefile <- st_read("02_transformed_data/EBSshelf/EBSshelf.shp")
crs_ref <- "+proj=longlat +ellps=WGS84  +datum=WGS84 +no_defs"
EBS <- st_transform(shapefile,crs=crs_ref)

p <- ggplot()+
  geom_sf(data = EBS,fill=NA, color="orange3",size=2) + 
  geom_sf(data = world,fill="black",color=NA) +
  coord_sf(xlim=c(-180,-155), ylim=c(54,63.5),expand = FALSE) 
p
ggsave("survey.png")


# -------------------------------------------------------------------------
devtools::install_github("MikkoVihtakari/ggOceanMapsData") # required by ggOceanMaps
devtools::install_github("MikkoVihtakari/ggOceanMaps")
library(ggOceanMaps)

dt <- data.frame(lon = c(-10, -10, 5,5), lat = c(40, 52, 52, 40))

q <- basemap(data = dt, bathymetry = TRUE,land.col = "black",land.border.col = NA, grid.col = NA) + labs(x="",y="",fill="") + theme(legend.position = "none")
q
ggsave("gulfGascogne.png")
