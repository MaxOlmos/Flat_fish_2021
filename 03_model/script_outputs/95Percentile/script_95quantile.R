
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
load("04_outputs/user_grid/EnvCxSeason/c(3,3,3,3,2,2)_grid1/fit.RData")
fit_f <- fit
dim(fit_f$Report$D_gct)
fit$data_list$n_i
fit_f$Report$Phi1_ik

# Calculate abundance because we want proportion
Spatial_List <- fit_f$spatial_list
a_k <- c(Spatial_List$a_xl[, 1], rep(0,Spatial_List$MeshList$anisotropic_mesh$n-nrow(Spatial_List$MeshList$loc_x)))
knot_i <- seq(1,length(a_k),1)
a_k <- as.data.frame(cbind(a_k,knot_i))
colnames(a_k) < c("a_k","knot")

# warm and cold seasons
cold <- c(6,7,8,9,10,11,12,13,17)
cold_year <- c(2006,2007,2008,2009,2010,2011,2012,2013,2017)

warm <- c(1,2,3,4,5,14,15,16,18)
warm_year <- c(2001,2002,2003,2004,2005,2014,2015,2016,2018)

n_knot <- fit_f$spatial_list$n_x 
name_effect <- c("ColdEarly","WarmEarly","ColdInt","WarmInt","ColdLate","WarmLate")
season <- c("Early","Int.","Late")
cold_season <- c(1,3,5)
warm_season <- c(2,4,6)


# -------------------------------------------------------------------------
# Abundances --------------------------------------------------------------
# -------------------------------------------------------------------------
dim(fit_f$Report$Phi1_sk)
D_cold_gut <- array(NA,dim=c(n_knot,length(cold_season),length(cold)))
dim(D_cold_gut)
D_warm_gut <- array(NA,dim=c(n_knot,length(warm_season),length(warm)))

for ( t in 1:length(cold)){
  for ( u in 1 :length(cold_season)){
    for ( g in 1:n_knot){
      D_cold_gut[g,u,t] <- fit_f$Report$D_gct[g,1,cold[t]] * a_k[g,1] * exp(fit_f$Report$Phi1_sk[g,cold_season[u]]) 
    }
  }
}

for ( t in 1:length(warm)){
  for ( u in 1 :length(warm_season)){
    for ( g in 1:n_knot){
      D_warm_gut[g,u,t] <- fit_f$Report$D_gct[g,1,warm[t]] *a_k[g,1] * exp(fit_f$Report$Phi1_sk[g,warm_season[u]])  
    }
  }
}
dim(fit_f$Report$D_gct)

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


rpoly_warm_new <- rpoly_warm %>% group_by(Season,EnvCond) %>%  summarize(geometry = st_union(geometry))
rpoly_cold_new <- rpoly_cold %>% group_by(Season,EnvCond) %>%  summarize(geometry = st_union(geometry))

rpoly_new <- bind_rows(rpoly_warm_new,rpoly_cold_new)

# define color
ramp <- colour_ramp(c("red", "blue"),alpha=TRUE)
ramp(seq(0, 1, length = 3))


#add survey
# -------------------------------------------------------------------------
load("04_outputs/EBS_grid/EnvCxSeason/survey/fit.RData")
fit_s <- fit
survey <- fit_s$data_frame %>% filter(t_i==2001)
concav_survey <-st_as_sf(survey, coords = c("Lon_i", "Lat_i"),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>% 
                concaveman::concaveman(concavity = 2.05)

# Plot  -------------------------------------------------------------------
p <- ggplot()+
  geom_point(survey,mapping=aes(x=Lon_i,y=Lat_i),size=0.8,color="lightgrey")+
  geom_sf(data = concav_survey,fill=NA, color="grey") + 
  geom_sf(data=rpoly_new,mapping=aes(color=EnvCond,fill=EnvCond),alpha=0.2,size=1.5)+theme_bw()+
  scale_color_manual(values=c("#1E90FFAA","#FF0000AA"))+
  scale_fill_manual(values=c( "#1E90FFAA","#FF0000AA"))+
  geom_sf(data = world,fill="black",color=NA) + 
  coord_sf(xlim=c(-180,-155), ylim=c(54,63.5),expand = FALSE) +
  facet_wrap( ~ Season) 

png(paste('Ab95_grid1','.png',sep=''), height =14 , width = 20, units = 'cm', res=500)#, res=500)
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
      D_cold_gut[g,u,t] <- exp(fit_f$Report$Phi1_sk[g,cold_season[u]]) 
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


rpoly_warm_new <- rpoly_warm %>% group_by(Season,EnvCond) %>%  summarize(geometry = st_union(geometry))
rpoly_cold_new <- rpoly_cold %>% group_by(Season,EnvCond) %>%  summarize(geometry = st_union(geometry))

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


# Load packages -----------------------------------------------------------
# -------------------------------------------------------------------------
library(VAST)
library(TMB)
library(dplyr)
library(rgdal)
library(raster)
library(pander)
library(splines)
library(FishStatsUtils)
library(devtools)
library(sf)


# meaning -----------------------------------------------------------------

#The 95% confidence interval means the probability that [pLo,pUp] contains the true cdf value is 0.95.
# cumulative distribution function (CDF) to find the percentage of samples with soil 
#test values at or below certain levels irrespective of their location.

# Fisheries ---------------------------------------------------------------
getwd()
# - fisheries CPUE --------------------------------------------------------
load("04_outputs/user_grid/EnvCxSeason/c(3,3,3,3,2,2)/fit.RData")
fit_f <- fit
dim(fit_f$Report$D_gct)


# warm and cold seasons
cold <- c(6,7,8,9,10,11,12,13,17)
cold_year <- c(2006,2007,2008,2009,2010,2011,2012,2013,2017)

warm <- c(1,2,3,4,5,14,15,16,18)
warm_year <- c(2001,2002,2003,2004,2005,2014,2015,2016,2018)

n_knot <- fit_f$spatial_list$n_x 
name_effect <- c("ColdEarly","WarmEarly","ColdInt","WarmInt","ColdLate","WarmLate")
season <- c("Early","Int.","Late")
cold_season <- c(1,3,5)
warm_season <- c(2,4,6)




# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Sex --------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
load("04_outputs/sex/fit_female.RData")
fit_female <- fit

load("04_outputs/sex/fit_male.RData")
fit_male <- fit

# Spatial resolution
Spatial_List <- fit_female$spatial_list
a_k <- c(Spatial_List$a_xl[, 1], rep(0,Spatial_List$MeshList$anisotropic_mesh$n-nrow(Spatial_List$MeshList$loc_x)))
knot_i <- seq(1,length(a_k),1)
a_k <- as.data.frame(cbind(a_k,knot_i))
colnames(a_k) < c("a_k","knot")

# rasterise ---------------------------------------------------------------
MapDetails_List <- FishStatsUtils:: make_map_info( "Region"=fit_female$settings$Region,spatial_list = fit_female$spatial_list, 
                                                   "NN_Extrap"=fit_female$spatial_list$PolygonList$NN_Extrap,                                                     
                                                   "Extrapolation_List"=fit_female$extrapolation_list ,fine_scale = FALSE)
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


# warm and cold seasons
# -------------------------------------------------------------------------

cold <- c(6,7,8,9,10,11,12,13,17)
cold_year <- c(2006,2007,2008,2009,2010,2011,2012,2013,2017)

warm <- c(1,2,3,4,5,14,15,16,18)
warm_year <- c(2001,2002,2003,2004,2005,2014,2015,2016,2018)

n_knot <- fit_female$spatial_list$n_x 
name_effect <- c("ColdEarly","WarmEarly","ColdInt","WarmInt","ColdLate","WarmLate")
season <- c("Early","Int.","Late")
cold_season <- c(1,3,5)
warm_season <- c(2,4,6)


D_cold_guts <- array(NA,dim=c(n_knot,length(cold_season),length(cold),2))
dim(D_cold_guts)
D_warm_guts <- array(NA,dim=c(n_knot,length(warm_season),length(warm),2))

# female =1
# male =2
for ( t in 1:length(cold)){
  for ( u in 1 :length(cold_season)){
    for ( g in 1:n_knot){
      D_cold_guts[g,u,t,1] <- fit_female$Report$D_gct[g,1,cold[t]] * a_k[g,1] * exp(fit_female$Report$Phi1_sk[g,cold_season[u]]) 
      D_cold_guts[g,u,t,2] <- fit_male$Report$D_gct[g,1,cold[t]] * a_k[g,1] * exp(fit_male$Report$Phi1_sk[g,cold_season[u]]) 
    }
  }
}

for ( t in 1:length(warm)){
  for ( u in 1 :length(warm_season)){
    for ( g in 1:n_knot){
      D_warm_guts[g,u,t,1] <- fit_female$Report$D_gct[g,1,warm[t]] *a_k[g,1] * exp(fit_female$Report$Phi1_sk[g,warm_season[u]])
      D_warm_guts[g,u,t,2] <- fit_male$Report$D_gct[g,1,warm[t]] *a_k[g,1] * exp(fit_male$Report$Phi1_sk[g,warm_season[u]])
    }
  }
}


# - Average accross all years ---------------------------------------------
D_cold_gu1 <- apply(D_cold_guts[,,,1],c(1,2),mean)
D_cold_gu2 <- apply(D_cold_guts[,,,2],c(1,2),mean)
D_warm_gu1 <- apply(D_warm_guts[,,,1],c(1,2),mean)
D_warm_gu2 <- apply(D_warm_guts[,,,2],c(1,2),mean)



# - 95% cumulative layers ---------------------------------------------------
# -------------------------------------------------------------------------
# cold
# - female
Y_g1 <- matrix(0,dim(D_cold_gu1)[1],dim(D_cold_gu1)[2])
order_g1 <- matrix(0,dim(D_cold_gu1)[1],dim(D_cold_gu1)[2])
cumsum_g1<- matrix(0,dim(D_cold_gu1)[1],dim(D_cold_gu1)[2])
threshold1<- matrix(0,dim(D_cold_gu1)[1],dim(D_cold_gu1)[2])
layer_t1<- matrix(0,dim(D_cold_gu1)[1],dim(D_cold_gu1)[2])
Y95_1<- matrix(0,dim(D_cold_gu1)[1],dim(D_cold_gu1)[2])

for (u in 1: ncol(D_cold_gu1))
{
  Y_g1[,u] = D_cold_gu1[,u]
  order_g1[,u] = order(Y_g1[,u],decreasing=TRUE)
  cumsum_g1[,u] = cumsum(Y_g1[order_g1[,u],u]) / sum(Y_g1[,u])
  threshold1[,u] = max( Y_g1[order_g1[cumsum_g1[,u]> 0.95],u] )
  layer_t1[,u] = ifelse( Y_g1[,u] > threshold1[,u], 1, 0 )
  Y95_1[,u] <- Y_g1[,u]*layer_t1[,u]
  
}
Y1_cold <- cbind(Y95_1,a_k[1:100,])
colnames(Y1_cold) =c("Early","Int","Late","a_k" ,"knot")


# - male
Y_g2 <- matrix(0,dim(D_cold_gu2)[1],dim(D_cold_gu2)[2])
order_g2 <- matrix(0,dim(D_cold_gu2)[1],dim(D_cold_gu2)[2])
cumsum_g2<- matrix(0,dim(D_cold_gu2)[1],dim(D_cold_gu2)[2])
threshold2<- matrix(0,dim(D_cold_gu2)[1],dim(D_cold_gu2)[2])
layer_t2<- matrix(0,dim(D_cold_gu2)[1],dim(D_cold_gu2)[2])
Y95_2<- matrix(0,dim(D_cold_gu2)[1],dim(D_cold_gu2)[2])

for (u in 1: ncol(D_cold_gu2))
{
  Y_g2[,u] = D_cold_gu2[,u]
  order_g2[,u] = order(Y_g2[,u],decreasing=TRUE)
  cumsum_g2[,u] = cumsum(Y_g2[order_g1[,u],u]) / sum(Y_g2[,u])
  threshold2[,u] = max( Y_g2[order_g2[cumsum_g2[,u]> 0.95],u] )
  layer_t2[,u] = ifelse( Y_g2[,u] > threshold2[,u], 1, 0 )
  Y95_2[,u] <- Y_g2[,u]*layer_t2[,u]
  
}


Y2_cold <- cbind(Y95_2,a_k[1:100,])
colnames(Y2_cold) =c("Early","Int","Late","a_k" ,"knot")


#warm
# - female
Y_g1 <- matrix(0,dim(D_warm_gu1)[1],dim(D_warm_gu1)[2])
order_g1 <- matrix(0,dim(D_warm_gu1)[1],dim(D_warm_gu1)[2])
cumsum_g1<- matrix(0,dim(D_warm_gu1)[1],dim(D_warm_gu1)[2])
threshold1<- matrix(0,dim(D_warm_gu1)[1],dim(D_warm_gu1)[2])
layer_t1<- matrix(0,dim(D_warm_gu1)[1],dim(D_warm_gu1)[2])
Y95_1<- matrix(0,dim(D_warm_gu1)[1],dim(D_warm_gu1)[2])

for (u in 1: ncol(D_warm_gu1))
{
  Y_g1[,u] = D_warm_gu1[,u]
  order_g1[,u] = order(Y_g1[,u],decreasing=TRUE)
  cumsum_g1[,u] = cumsum(Y_g1[order_g1[,u],u]) / sum(Y_g1[,u])
  threshold1[,u] = max( Y_g1[order_g1[cumsum_g1[,u]> 0.95],u] )
  layer_t1[,u] = ifelse( Y_g1[,u] > threshold1[,u], 1, 0 )
  Y95_1[,u] <- Y_g1[,u]*layer_t1[,u]
  
}
Y1_warm <- cbind(Y95_1,a_k[1:100,])
colnames(Y1_warm) =c("Early","Int","Late","a_k" ,"knot")


# - male
Y_g2 <- matrix(0,dim(D_warm_gu2)[1],dim(D_warm_gu2)[2])
order_g2 <- matrix(0,dim(D_warm_gu2)[1],dim(D_warm_gu2)[2])
cumsum_g2<- matrix(0,dim(D_warm_gu2)[1],dim(D_warm_gu2)[2])
threshold2<- matrix(0,dim(D_warm_gu2)[1],dim(D_warm_gu2)[2])
layer_t2<- matrix(0,dim(D_warm_gu2)[1],dim(D_warm_gu2)[2])
Y95_2<- matrix(0,dim(D_warm_gu2)[1],dim(D_warm_gu2)[2])

for (u in 1: ncol(D_warm_gu2))
{
  Y_g2[,u] = D_warm_gu2[,u]
  order_g2[,u] = order(Y_g2[,u],decreasing=TRUE)
  cumsum_g2[,u] = cumsum(Y_g2[order_g1[,u],u]) / sum(Y_g2[,u])
  threshold2[,u] = max( Y_g2[order_g2[cumsum_g2[,u]> 0.95],u] )
  layer_t2[,u] = ifelse( Y_g2[,u] > threshold2[,u], 1, 0 )
  Y95_2[,u] <- Y_g2[,u]*layer_t2[,u]
  
}


Y2_warm <- cbind(Y95_2,a_k[1:100,])
colnames(Y2_warm) =c("Early","Int","Late","a_k" ,"knot")


##############################################################
# early
Poly <- function(Y,Cond, sex){

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

rpoly <- bind_rows(r_poly1,r_poly2,r_poly3) %>% mutate(EnvCond=Cond, Sex = sex)
return(rpoly)

}

rpoly_female_warm <- Poly(Y1_warm,"warm","female") 
rpoly_female_cold <-  Poly(Y1_cold,"cold","female")  
rpoly_male_warm <- Poly(Y2_warm,"warm","male") 
rpoly_male_cold <- Poly(Y2_cold,"cold","male") 



# Maps --------------------------------------------------------------------
# -------------------------------------------------------------------------
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)
library(ggplot2)
library(concaveman)

world <- ne_countries(scale = "medium", returnclass = "sf")

rpoly_female_warm_new <- rpoly_female_warm %>% group_by(Season,EnvCond,Sex) %>%  summarize(geometry = st_union(geometry))
rpoly_female_cold_new <-  rpoly_female_cold %>% group_by(Season,EnvCond,Sex) %>%  summarize(geometry = st_union(geometry)) 
rpoly_male_warm_new <- rpoly_male_warm %>% group_by(Season,EnvCond,Sex) %>%  summarize(geometry = st_union(geometry))
rpoly_male_cold_new <- rpoly_male_cold %>% group_by(Season,EnvCond,Sex) %>%  summarize(geometry = st_union(geometry))


rpoly_new <- bind_rows(rpoly_female_warm_new,rpoly_female_cold_new,rpoly_male_warm_new,rpoly_male_cold_new)

# define color
ramp <- colour_ramp(c("red", "blue"),alpha=TRUE)
ramp(seq(0, 1, length = 3))


#add survey
# -------------------------------------------------------------------------
load("04_outputs/EBS_grid/EnvCxSeason/survey/fit.RData")
fit_s <- fit
survey <- fit_s$data_frame %>% filter(t_i==2001)
concav_survey <-st_as_sf(survey, coords = c("Lon_i", "Lat_i"),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>% 
  concaveman::concaveman(concavity = 2.05)

# Plot  -------------------------------------------------------------------
#EnvCond_vs_sex

p <- ggplot()+
  geom_point(survey,mapping=aes(x=Lon_i,y=Lat_i),size=0.8,color="lightgrey")+
  geom_sf(data = concav_survey,fill=NA, color="grey") + 
  geom_sf(data=rpoly_new,mapping=aes(color=EnvCond,fill=EnvCond),alpha=0.2,size=1.5)+theme_bw()+
  scale_color_manual(values=c("#1E90FFAA","#FF0000AA"))+
  scale_fill_manual(values=c( "#1E90FFAA","#FF0000AA"))+
  geom_sf(data = world,fill="black",color=NA) + 
  coord_sf(xlim=c(-180,-155), ylim=c(54,63.5),expand = FALSE) +
  facet_grid( Sex ~ Season) 
p
png(paste('EnvCond_vs_sex','.png',sep=''), height =14 , width = 20, units = 'cm', res=500)#, res=500)
p
dev.off()

#sex_vs_EnvCond
p <- ggplot()+
  geom_point(survey,mapping=aes(x=Lon_i,y=Lat_i),size=0.8,color="lightgrey")+
  geom_sf(data = concav_survey,fill=NA, color="grey") + 
  geom_sf(data=rpoly_new,mapping=aes(color=Sex,fill=Sex),alpha=0.2,size=1.5)+theme_bw()+
  scale_color_manual(values=c("#009900AA","#FF9933AA"))+
  scale_fill_manual(values=c( "#009900AA","#FF9933AA"))+
  geom_sf(data = world,fill="black",color=NA) + 
  coord_sf(xlim=c(-180,-155), ylim=c(54,63.5),expand = FALSE) +
  facet_grid( EnvCond ~ Season) 
p
png(paste('sex_vs_EnvCond','.png',sep=''), height =14 , width = 20, units = 'cm', res=500)#, res=500)
p
dev.off()



