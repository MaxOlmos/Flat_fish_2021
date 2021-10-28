library(metR)
library(ggplot2)
library(VAST)                           # 3.8.0
library(ggplot2)                        # 2.10.0
library(dplyr)
library(tidyr)
theme_set(theme_bw())
library(concaveman)
library(sf)
dev.off()
getwd()


# Load bathy --------------------------------------------------------------
bathy.dat<-read.table('C:/Users/Maxime/Documents/Git/Flat_fish_2021/01_data/bathyEBS/BeringDepth.txt',sep='') 
names(bathy.dat)<-c('lon','lat','depth')
head(bathy.dat)

range(fit$data_frame$Lat_i)
range(fit$data_frame$Lon_i)


# Load spatial limits -----------------------------------------------------

plot(c(fit$data_frame$Lon_i), fit$data_frame$Lat_i)

LL <- locator()
LL2 <- cbind(LL$x, LL$y)
LL2 <- as.data.frame(LL2)
load("02_transformed_data/bathyEBS/LL2.RData")
bathy.dat <- bathy.dat %>% filter(lon < -155, lon > -179,lat<61, lat>54)
polygon_bathy <- st_as_sf(as.data.frame(LL2),coords = c("V1", "V2"))%>% 
  concaveman::concaveman(concavity = 5/4)
point <- st_as_sf(bathy.dat,coords = c("lon","lat"))

bathy.dat2 <- st_intersection(point,polygon_bathy$polygons )  
bathy.dat3 <- cbind(st_coordinates(bathy.dat2),bathy.dat2$depth)
bathy.dat3 <- as.data.frame(bathy.dat3)
colnames(bathy.dat3) <- c("X","Y","depth")



# prepare facet ------------------------------------------------------------
#try with fine scale = TRUE
## Remake map list locally for recreating plots
mdl <- make_map_info(Region ="user",
                     spatial_list = fit$spatial_list,
                     Extrapolation_List = fit$extrapolation_list)

mdl$Xlim <- c(mdl$Xlim[1],mdl$Xlim[2]+5)
mdl$Ylim <- c(mdl$Ylim[1],mdl$Ylim[2]+1)

years <- unique(fit$data_frame$t_i)
nyrs <- length(years)

name_effect <- c("ColdEarly","ColdInt","ColdLate","WarmEarly","WarmInt","WarmLate")
neffect <- length(name_effect)

## quick dirty AK map
ak_map <- subset(map_data("world"), region=='USA' & subregion=='Alaska')

## Have to duplicate it for each year so can facet below
ak_map_effect <- cbind(ak_map[rep(1:nrow(ak_map), times=neffect),],
                Effect=rep(name_effect, each=nrow(ak_map)))

ak_map <- cbind(ak_map[rep(1:nrow(ak_map), times=nyrs),],
                Year=rep(years, each=nrow(ak_map)))
gmap_effect <- ggplot(ak_map_effect, aes(x = long, y = lat, group = group)) +
   geom_polygon(fill="black", colour = "white") +
  scale_color_viridis_c(option = "viridis") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.spacing.x=unit(0, "lines"),
        panel.spacing.y=unit(0, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ) +
  coord_cartesian(xlim=mdl$Xlim, ylim=mdl$Ylim)


names(fit$Report)[grepl('Phi', x=names(fit$Report))]

phi_sk <-fit$Report$Phi1_sk # drop the category

name_effect <- c("ColdEarly","WarmEarly","ColdInt","WarmInt","ColdLate","WarmLate")

dimnames(phi_sk) <- list(cell=1:nrow(phi_sk), year=name_effect)
## tidy way of doing this, reshape2::melt() does
## it cleanly but is deprecated
phi_sk <- phi_sk %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "cell") %>%
  pivot_longer(-cell, names_to = "Effect", values_to='phi')
phi<- merge(phi_sk, mdl$PlotDF, by.x='cell', by.y='x2i')


phi$Effect <- factor(phi$Effect, levels = c("ColdEarly","ColdInt","ColdLate","WarmEarly","WarmInt","WarmLate"))
g <- gmap_effect+
  theme(panel.grid = element_line(color = "red"),
        panel.ontop = TRUE)+
  geom_point(data=phi, aes(Lon, Lat, color=(phi), group=NULL),
             ## These settings are necessary to avoid
             ## overlplotting which is a problem here. May need
             ## to be tweaked further.
             size=3,
             stroke=0,shape=16) + facet_wrap('Effect')  +
  theme_bw()+ geom_contour( data=bathy.dat,aes(x = lon, y = lat,z = depth,group=NULL),breaks  = c(-200,-100,-50),col="green4",size=1)

g






ggplot(ak_map_effect, aes(x = long, y = lat, group = group)) +

  geom_point(data=phi, aes(Lon, Lat, color=phi, group=NULL),
             ## These settings are necessary to avoid
             ## overlplotting which is a problem here. May need
             ## to be tweaked further.
             size=4,
             stroke=0,shape=16) + facet_wrap('Effect')  +
  theme_bw()+ geom_contour( data=bathy.dat3,aes(x = X, y = Y,z = depth,group=NULL),breaks  = c(-200,-100,-50),col="black",size=1)+
  
  geom_polygon(fill="black", colour = "white") +
  scale_color_viridis_c(option = "magma") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.spacing.x=unit(0, "lines"),
        panel.spacing.y=unit(0, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ) +
  coord_cartesian(xlim=mdl$Xlim, ylim=mdl$Ylim)





names(fit$Report)[grepl('_gc|_gct', x=names(fit$Report))]
dim(fit$Report$D_gct[,1,])
D_gt <- fit$Report$D_gct[,1,] # drop the category
dimnames(D_gt) <- list(cell=1:nrow(D_gt), year=years)
## tidy way of doing this, reshape2::melt() does
## it cleanly but is deprecated
D_gt <- D_gt %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "cell") %>%
  pivot_longer(-cell, names_to = "Year", values_to='D')

unique(D_gt$cell)
D <- merge(D_gt, mdl$PlotDF, by.x='cell', by.y='x2i')
unique(D$cell)

g <- gmap +
  geom_point(data=D, aes(Lon, Lat, color=log(D), group=NULL),
             ## These settings are necessary to avoid
             ## overlplotting which is a problem here. May need
             ## to be tweaked further.
             
             size=2.5, stroke=0,shape=16) + facet_wrap('Year')+ 
  theme_bw()
g

## Or you can do it in base R, building your own palette and
## looping through years as needed.