rm(list=ls())

# Load packages -----------------------------------------------------------
# -------------------------------------------------------------------------
library(VAST)
library(TMB)
library(dplyr)
library(rgdal)
library(raster)
library(pander)
library(splines)
library(devtools)
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")
library("scales")
library(ggplot2)
library(akgfmaps)

# - 95% cumulative layers ---------------------------------------------------
# -------------------------------------------------------------------------
# for example here I am calculating the 95% percentiles for the average density accross cold years : D_cold_gu
# where g is for the spatial resolution (cells) and u is for the different seasons (early, intermediate, late)
load("D_cold_gu.Rdata")
head(D_cold_gu)

Y_g <- matrix(0,dim(D_cold_gu)[1],dim(D_cold_gu)[2])
order_g <- matrix(0,dim(D_cold_gu)[1],dim(D_cold_gu)[2])
cumsum_g<- matrix(0,dim(D_cold_gu)[1],dim(D_cold_gu)[2])
threshold<- matrix(0,dim(D_cold_gu)[1],dim(D_cold_gu)[2])
layer_t<- matrix(0,dim(D_cold_gu)[1],dim(D_cold_gu)[2])
Y95<- matrix(0,dim(D_cold_gu)[1],dim(D_cold_gu)[2])

for (u in 1: ncol(D_cold_gu))
{
  Y_g[,u] = D_cold_gu[,u]
  order_g[,u] = order(Y_g[,u],decreasing=TRUE)
  cumsum_g[,u] = cumsum(Y_g[order_g[,u],u]) / sum(Y_g[,u])
  threshold[,u] = max( Y_g[order_g[cumsum_g[,u]> 0.95],u] )
  layer_t[,u] = ifelse( Y_g[,u] > threshold[,u], 1, 0 )
  Y95[,u] <- Y_g[,u]*layer_t[,u]
  
}

Y <-Y95
colnames(Y) =c("Early","Int","Late")

# Extrapolation grid : EG -------------------------------------------------
# Rasterise ---------------------------------------------------------------

load( file="EG.Rdata") # here is for example but use yours (if needed I can show you how to build a extrapolation grid)
plot(EG)

PlotDF1 <- as.data.frame(EG)
coords2 = cbind(PlotDF1$Lon,PlotDF1$Lat)
coords2 = coords2[coords2[,1]<0,]
P2 = SpatialPoints(coords2)

rast <- raster(ncol=50,nrow=40) # feel free to change the dimension of the raster
extent(rast) <- extent(coords2)

# Rasterise 95% percential of the spatial distribution for each ---------
# season ----------------------------------------------------------------

# early
P2$data =  Y[,1]
rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
r_poly1 <- rasterToPolygons(rast.temp, function(x ){x >0}, dissolve = TRUE) %>% 
  st_as_sf() %>% mutate(Season="Early", value=layer) 

# int
P2$data =  Y[,2]
rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
r_poly2 <- rasterToPolygons(rast.temp, function(x ){x >0}, dissolve = TRUE) %>% 
  st_as_sf() %>% mutate(Season="Int", value=layer) 

# late
P2$data =  Y[,3]
rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
r_poly3 <- rasterToPolygons(rast.temp, function(x ){x >0}, dissolve = TRUE) %>% 
  st_as_sf()%>% mutate(Season="Late", value=layer) 

rpoly_cold <- bind_rows(r_poly1,r_poly2,r_poly3) %>% mutate(EnvCond="cold")

# Define a polygon (an area) which include 95% percentile
rpoly_cold_merged<- rpoly_cold %>% dplyr::group_by(Season,EnvCond) %>%  dplyr::summarize(geometry = st_union(geometry))


# Maps --------------------------------------------------------------------
# -------------------------------------------------------------------------

#- define world
world <- ne_countries(scale = "medium", returnclass = "sf")

#- define color
ramp <- colour_ramp(c("red", "blue"),alpha=TRUE)
ramp(seq(0, 1, length = 3))

#- bathy
crs_ref <- "+proj=longlat +ellps=WGS84  +datum=WGS84 +no_defs"
sebs_layers <- akgfmaps::get_base_layers(select.region = "sebs",
                                         set.crs =crs_ref)

# Plot  -------------------------------------------------------------------
# - 1 Rasterwith values colored
# here tje plot by itself does not make any sense in term of ecology but it is just to show you how to rasterise (plot2 makes more sense))

plot <- ggplot()+
  geom_sf(data=rpoly_cold, aes(fill=log(value)),col=NA)+
  scale_fill_viridis() +
  geom_sf(data = sebs_layers$bathymetry,size=0.2) + 
  geom_sf(data = world,fill="black",color=NA) + 
  coord_sf(xlim=c(-180,-155), ylim=c(54,63.5),expand = FALSE)+
  theme_bw()

plot1 <- plot +  facet_wrap( ~ Season) 

ggsave('plot1.png',plot1)


# - 2 Raster with only 95% areas 
plot <- ggplot()+
  geom_sf(data = sebs_layers$bathymetry,size=0.2) + 
  geom_sf(data=rpoly_cold_merged,mapping=aes(color=EnvCond,fill=EnvCond),alpha=0.2,size=1.5)+
  scale_color_manual(values=c("#1E90FFAA"))+
  scale_fill_manual(values=c( "#1E90FFAA"))+
  geom_sf(data = world,fill="black",color=NA) + 
  coord_sf(xlim=c(-180,-155), ylim=c(54,63.5),expand = FALSE) +
  theme_bw()+
  facet_wrap( ~ Season)

plot2 <- plot +  facet_wrap( ~ Season) 

ggsave('plot2.png',plot2)

