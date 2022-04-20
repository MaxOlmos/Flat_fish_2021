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
#library(FishStatsUtils)
library(devtools)
library(sf)


# meaning -----------------------------------------------------------------

#The 95% confidence interval means the probability that [pLo,pUp] contains the true cdf value is 0.95.
# F function (CDF) to find the percentage of samples with soil 
#test values at or below certain levels irrespective of their location.

# Fisheries ---------------------------------------------------------------
getwd()
# - fisheries CPUE --------------------------------------------------------
work_dir <- "04_outputs/user_grid_outEBS/no_sex/finescale/c(3,3,3,3,2,2)/no2000/av_season/"
load(paste0(work_dir,"fit.RData"))
fit_f <- fit
dim(fit_f$Report$D_gct)

# warm and cold seasons
#cold <- c(1,7,8,9,10,11,12,13,14,18)
#cold_year <- c(2000,2006,2007,2008,2009,2010,2011,2012,2013,2017)

#warm <- c(2,3,4,5,6,15,16,17,19,20)
#warm_year <- c(2001,2002,2003,2004,2005,2014,2015,2016,2018,2019)

cold <- c(6,7,8,9,10,11,12,13,17)
cold_year <- c(2006,2007,2008,2009,2010,2011,2012,2013,2017)

warm <- c(1,2,3,4,5,14,15,16,18,19)
warm_year <- c(2001,2002,2003,2004,2005,2014,2015,2016,2018,2019)


name_effect <- c("ColdEarly","WarmEarly","ColdInt","WarmInt","ColdLate","WarmLate")
season <- c("Early","Int.","Late")
cold_season <- c(1,3,5)
warm_season <- c(2,4,6)

sd1 <- summary(fit$parameter_estimates$SD)[,1]
sd2 <- fit$parameter_estimates$SD[,2]
as.vector(sd1[which(names(sd1)=="Phi1_sk")])
unique(names(sd1))


# -------------------------------------------------------------------------
# Abundances --------------------------------------------------------------
# -------------------------------------------------------------------------


D_cold_gut <- array(NA,dim=c(dim(fit_f$Report$Phi1_gk)[1],length(cold_season),length(cold)))
D_warm_gut <- array(NA,dim=c(dim(fit_f$Report$Phi1_gk)[1],length(warm_season),length(warm)))

for ( t in 1:length(cold)){
  for ( u in 1 :length(cold_season)){
    for ( g in 1:dim(fit_f$Report$Phi1_gk)[1]){
      D_cold_gut[g,u,t] <- fit_f$Report$D_gct[g,,cold[t]] * exp(fit_f$Report$Phi1_gk[g,cold_season[u]]) 
    }
  }
}

for ( t in 1:length(warm)){
  for ( u in 1 :length(warm_season)){
    for ( g in 1:dim(fit_f$Report$Phi1_gk)[1]){
      D_warm_gut[g,u,t] <- fit_f$Report$D_gct[g,1,warm[t]] * exp(fit_f$Report$Phi1_gk[g,warm_season[u]])  
    }
  }
}


# - Average accross all years ---------------------------------------------
dim(D_cold_gut)
D_cold_gu <- apply(D_cold_gut ,c(1,2),mean)
D_warm_gu <- apply(D_warm_gut,c(1,2),mean)



# - 95% cumulative layers ---------------------------------------------------
# -------------------------------------------------------------------------
#cold

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

#warm

Y_g <- matrix(0,dim(D_warm_gu)[1],dim(D_warm_gu)[2])
order_g <- matrix(0,dim(D_warm_gu)[1],dim(D_warm_gu)[2])
cumsum_g<- matrix(0,dim(D_warm_gu)[1],dim(D_warm_gu)[2])
threshold<- matrix(0,dim(D_warm_gu)[1],dim(D_warm_gu)[2])
layer_t<- matrix(0,dim(D_warm_gu)[1],dim(D_warm_gu)[2])

for (u in 1: ncol(D_warm_gu))
{
  
  Y_g[,u] = D_warm_gu[,u]
  order_g[,u] = order(Y_g[,u],decreasing=TRUE)
  cumsum_g[,u] = cumsum(Y_g[order_g[,u],u]) / sum(Y_g[,u])
  threshold[,u] = max( Y_g[order_g[cumsum_g[,u]> 0.95],u] )
  layer_t[,u] = ifelse( Y_g[,u] > threshold[,u], 1, 0 )
  Y95[,u] <- Y_g[,u]*layer_t[,u]
  
}

Y <- Y95
colnames(Y) =c("Early","Int","Late")

# rasterise ---------------------------------------------------------------
MapDetails_List <- FishStatsUtils:: make_map_info( "Region"=fit_f$settings$Region,spatial_list = fit_f$spatial_list, "NN_Extrap"=fit_f$spatial_list$PolygonList$NN_Extrap,                                                     
                                                   "Extrapolation_List"=fit_f$extrapolation_list ,fine_scale = FALSE)
x <- as.data.frame(MapDetails_List$PlotDF)
PlotDF1 = subset(x,Include==TRUE)

unique(x$Include)
PlotDF1 <- as.data.frame(PlotDF1)
coords2 = cbind(PlotDF1$Lon,PlotDF1$Lat)
coords2 = coords2[coords2[,1]<0,]
P2 = SpatialPoints(coords2)
P3 = SpatialPoints(coords2)

rast <- raster(ncol=50,nrow=40)
extent(rast) <- extent(coords2)

##############################################################

# early
P2$data =  Y[PlotDF1$x2i,1]
rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
r_poly1 <- rasterToPolygons(rast.temp, function(x ){x >0}, dissolve = TRUE) %>% 
  st_as_sf() %>% mutate(Season="Early") 

# int
P2$data =  Y[PlotDF1$x2i,2]
rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
r_poly2 <- rasterToPolygons(rast.temp, function(x ){x >0}, dissolve = TRUE) %>% 
  st_as_sf() %>% mutate(Season="Int") 

# late
P2$data =  Y[PlotDF1$x2i,3]
rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
r_poly3 <- rasterToPolygons(rast.temp, function(x ){x >0}, dissolve = TRUE) %>% 
  st_as_sf()%>% mutate(Season="Late") 


rpoly_warm <- bind_rows(r_poly1,r_poly2,r_poly3) %>% mutate(EnvCond="warm")
rpoly_cold <- bind_rows(r_poly1,r_poly2,r_poly3) %>% mutate(EnvCond="cold")



# Maps --------------------------------------------------------------------
# -------------------------------------------------------------------------
library("rnaturalearth")
library("rnaturalearthdata")
library("scales")
library(ggplot2)


world <- ne_countries(scale = "medium", returnclass = "sf")
#ggplot()+geom_sf(data=rpoly,mapping=aes(fill=Season),color = NA,alpha=0.3 )+theme_bw()+
#geom_sf(data = world,fill="black") + 
# coord_sf(xlim=c(-180,-155), ylim=c(54,63.5), expand = FALSE)


rpoly_warm_new <- rpoly_warm %>% dplyr::group_by(Season,EnvCond) %>%  dplyr::summarize(geometry = st_union(geometry))
rpoly_cold_new <- rpoly_cold %>% dplyr::group_by(Season,EnvCond) %>%  dplyr::summarize(geometry = st_union(geometry))

rpoly_new <- bind_rows(rpoly_warm_new,rpoly_cold_new)

# define color
ramp <- colour_ramp(c("red", "blue"),alpha=TRUE)
ramp(seq(0, 1, length = 3))


#add survey
# -------------------------------------------------------------------------
shapefile <- st_read("02_transformed_data/EBSshelf/EBSshelf.shp")
crs_ref <- "+proj=longlat +ellps=WGS84  +datum=WGS84 +no_defs"
EBS <- st_transform(shapefile,crs=crs_ref)

## bathy
library(akgfmaps)
sebs_layers <- akgfmaps::get_base_layers(select.region = "sebs",
                                         set.crs =crs_ref)

# Plot  -------------------------------------------------------------------
p <- ggplot()+
  #geom_point(survey,mapping=aes(x=Lon_i,y=Lat_i),size=0.8,color="lightgrey")+
  #geom_sf(data = concav_survey,fill=NA, color="grey") + 
  geom_sf(data=EBS,fill=NA, color="brown")+
  geom_sf(data = sebs_layers$bathymetry,size=0.2) + 
  geom_sf(data=rpoly_new,mapping=aes(color=EnvCond,fill=EnvCond),alpha=0.2,size=1.5)+theme_bw()+
  scale_color_manual(values=c("#1E90FFAA","#FF0000AA"))+
  scale_fill_manual(values=c( "#1E90FFAA","#FF0000AA"))+
  geom_sf(data = world,fill="black",color=NA) + 
  coord_sf(xlim=c(-180,-155), ylim=c(54,63.5),expand = FALSE) +
  facet_wrap( ~ Season) 

png(paste0(work_dir,'Ab95_grid1','.png',sep=''), height =14 , width = 20, units = 'cm', res=500)#, res=500)
p
dev.off()


# -------------------------------------------------------------------------
# ------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Spatial effect ----------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

D_cold_gut <- array(NA,dim=c(n_knot,length(cold_season),length(cold)))
dim(D_cold_gut)
D_warm_gut <- array(NA,dim=c(n_knot,length(warm_season),length(warm)))

for ( t in 1:length(cold)){
  for ( u in 1 :length(cold_season)){
    for ( g in 1:n_knot){
      D_cold_gut[g,u,t] <- exp(fit_f$Report$Phi1_gk[g,cold_season[u]]) 
    }
  }
}

for ( t in 1:length(warm)){
  for ( u in 1 :length(warm_season)){
    for ( g in 1:n_knot){
      D_warm_gut[g,u,t] <-  exp(fit_f$Report$Phi1_sk[g,warm_season[u]])  
    }
  }
}


# - Average accross all years ---------------------------------------------
D_cold_gu <- apply(D_cold_gut,c(1,2),mean)
D_warm_gu <- apply(D_warm_gut,c(1,2),mean)



# - 95% cumulative layers ---------------------------------------------------
# -------------------------------------------------------------------------
#cold

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


Y <- cbind(Y95,a_k[1:100,])
colnames(Y) =c("Early","Int","Late","a_k" ,"knot")

#warm

Y_g <- matrix(0,dim(D_warm_gu)[1],dim(D_warm_gu)[2])
order_g <- matrix(0,dim(D_warm_gu)[1],dim(D_warm_gu)[2])
cumsum_g<- matrix(0,dim(D_warm_gu)[1],dim(D_warm_gu)[2])
threshold<- matrix(0,dim(D_warm_gu)[1],dim(D_warm_gu)[2])
layer_t<- matrix(0,dim(D_warm_gu)[1],dim(D_warm_gu)[2])

for (u in 1: ncol(D_warm_gu))
{
  
  Y_g[,u] = D_warm_gu[,u]
  order_g[,u] = order(Y_g[,u],decreasing=TRUE)
  cumsum_g[,u] = cumsum(Y_g[order_g[,u],u]) / sum(Y_g[,u])
  threshold[,u] = max( Y_g[order_g[cumsum_g[,u]> 0.95],u] )
  layer_t[,u] = ifelse( Y_g[,u] > threshold[,u], 1, 0 )
  Y95[,u] <- Y_g[,u]*layer_t[,u]
  
}

Y <- cbind(Y95,a_k[1:100,])
colnames(Y) =c("Early","Int","Late","a_k" ,"knot")

# rasterise ---------------------------------------------------------------
MapDetails_List <- FishStatsUtils:: make_map_info( "Region"=fit_f$settings$Region,spatial_list = fit_f$spatial_list, "NN_Extrap"=fit_f$spatial_list$PolygonList$NN_Extrap,                                                     
                                                   "Extrapolation_List"=fit_f$extrapolation_list ,fine_scale = FALSE)
x <- as.data.frame(MapDetails_List$PlotDF)
PlotDF1 = subset(x,Include==TRUE)

unique(x$Include)
PlotDF1 <- as.data.frame(PlotDF1)
coords2 = cbind(PlotDF1$Lon,PlotDF1$Lat)
coords2 = coords2[coords2[,1]<0,]
P2 = SpatialPoints(coords2)
P3 = SpatialPoints(coords2)

rast <- raster(ncol=50,nrow=40)
extent(rast) <- extent(coords2)

##############################################################

# early
P2$data =  Y[PlotDF1$x2i,1]
rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
r_poly1 <- rasterToPolygons(rast.temp, function(x ){x >0}, dissolve = TRUE) %>% 
  st_as_sf() %>% mutate(Season="Early") 

# int
P2$data =  Y[PlotDF1$x2i,2]
rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
r_poly2 <- rasterToPolygons(rast.temp, function(x ){x >0}, dissolve = TRUE) %>% 
  st_as_sf() %>% mutate(Season="Int") 

# late
P2$data =  Y[PlotDF1$x2i,3]
rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
r_poly3 <- rasterToPolygons(rast.temp, function(x ){x >0}, dissolve = TRUE) %>% 
  st_as_sf()%>% mutate(Season="Late") 


rpoly_warm <- bind_rows(r_poly1,r_poly2,r_poly3) %>% mutate(EnvCond="warm")
rpoly_cold <- bind_rows(r_poly1,r_poly2,r_poly3) %>% mutate(EnvCond="cold")



# Maps --------------------------------------------------------------------
# -------------------------------------------------------------------------
library("rnaturalearth")
library("rnaturalearthdata")
library("scales")
library(ggplot2)


world <- ne_countries(scale = "medium", returnclass = "sf")
#ggplot()+geom_sf(data=rpoly,mapping=aes(fill=Season),color = NA,alpha=0.3 )+theme_bw()+
#geom_sf(data = world,fill="black") + 
# coord_sf(xlim=c(-180,-155), ylim=c(54,63.5), expand = FALSE)


rpoly_warm_new <- rpoly_warm %>% dplyr::group_by(Season,EnvCond) %>% dplyr:: summarize(geometry = st_union(geometry))
rpoly_cold_new <- rpoly_cold %>% dplyr::group_by(Season,EnvCond) %>% dplyr:: summarize(geometry = st_union(geometry))

rpoly_new <- bind_rows(rpoly_warm_new,rpoly_cold_new)

# define color
ramp <- colour_ramp(c("red", "blue"),alpha=TRUE)
ramp(seq(0, 1, length = 3))


#add survey
# -------------------------------------------------------------------------
#load("04_outputs/EBS_grid/EnvCxSeason/survey/fit.RData")
#fit_s <- fit
#survey <- fit_s$data_frame %>% filter(t_i==2001)
#concav_survey <-st_as_sf(survey, coords = c("Lon_i", "Lat_i"),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>% 
# concaveman::concaveman(concavity = 2.05)

# Plot  -------------------------------------------------------------------
p <- ggplot()+
  # geom_point(survey,mapping=aes(x=Lon_i,y=Lat_i),size=0.8,color="lightgrey")+
  #  geom_sf(data = concav_survey,fill=NA, color="grey") + 
  geom_sf(data=rpoly_new,mapping=aes(color=EnvCond,fill=EnvCond),alpha=0.2,size=1.5)+theme_bw()+
  scale_color_manual(values=c("#1E90FFAA","#FF0000AA"))+
  scale_fill_manual(values=c( "#1E90FFAA","#FF0000AA"))+
  geom_sf(data = world,fill="black",color=NA) + 
  coord_sf(xlim=c(-180,-155), ylim=c(54,63.5),expand = FALSE) +
  facet_wrap( ~ Season) 
p
png(paste('Effect95','.png',sep=''), height =14 , width = 20, units = 'cm', res=500)#, res=500)
p
dev.off()


# -------------------------------------------------------------------------
library(coldpool)
load("C:/Users/Maxime/Documents/R/win-library/4.1/coldpool/R/sysdata.rda")

?coldpool:::cold_pool_index

# View loaded cold pool index data frame
cold_pool_index
ggplot() +geom_line(data=cold_pool_index, mapping=aes(x=YEAR,y=MEAN_GEAR_TEMPERATURE))

cpi_df

cpi_df <- cold_pool_index %>% 
  dplyr::filter(YEAR <= 2020, YEAR > 1999) %>%
  dplyr::mutate(group = YEAR < 2021,
                var = "Cold Pool Index",
                cpi_rank = rank(AREA_LTE2_KM2))

cpi_sd <- sd(cpi_df$AREA_LTE2_KM2/1000)
cpi_mean <- mean(cpi_df$AREA_LTE2_KM2/1000)
coeff=0.01
plot_cpi_timeseries <- ggplot(data = cpi_df %>% filter(YEAR>1988,YEAR<2021),
                              mapping = aes(x = YEAR) )+
  # geom_point(size = rel(2),
  #           color = "#000040") +
  #geom_line( aes(y = MEAN_GEAR_TEMPERATURE,
  #          group = group
  #   ), color = "blue",size=1.3) +
  geom_line( aes(y = AREA_LTE2_KM2/1000,
                 group = group,
                 label = cpi_rank), color = "blue",size=1.3) +
  geom_hline(yintercept = cpi_mean,
             linetype = 3,
             color = "blue") +
  geom_hline(yintercept = cpi_mean + c(cpi_sd, -1*cpi_sd),
             linetype = 4,
             color = "blue") +
  
  scale_x_continuous(name = "Year",
                     breaks = seq(1980,2020,4)) +
  #facet_wrap(~var) +
  theme_bw()




