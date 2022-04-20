# Load packages -----------------------------------------------------------
# -------------------------------------------------------------------------
library(VAST)
library(TMB)
library(tidyr)
library(sf)
library(rgdal)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library("rnaturalearth")
library("rnaturalearthdata")
library(reshape2)

getwd()

# Fisheries ---------------------------------------------------------------
load("04_outputs/user_grid/EnvCxSeason/c(3,3,3,3,2,2)_grid1/fit.RData")
fit_f <- fit


# Overlap steps -----------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# - 1 Get the shapefile for the EBS survey e.g. here:  --------------------
#https://github.com/James-Thorson-NOAA/FishStatsUtils/tree/main/inst/region_shapefiles/EBSshelf
# -------------------------------------------------------------------------

world <- ne_countries(scale = "medium", returnclass = "sf")

shapefile <- st_read(paste0(getwd(),"/02_transformed_data/EBSshelf/EBSshelf.shp"))
plot(shapefile$geometry)
st_crs(shapefile)

shapefile <- st_transform(shapefile, "+proj=longlat +datum=WGS84")
# plot
p <- ggplot() + 
  geom_sf(data = shapefile) + 
  ggtitle("EBS") 
p


# Survey points -----------------------------------------------------------
# -------------------------------------------------------------------------

load("04_outputs/EBS_grid/EnvCxSeason/survey/fit.RData")
fit_s <- fit
survey <- fit_s$data_frame 


# -------------------------------------------------------------------------
# - 2 Get the extrapolation-grid you’re using, fit$spatial_list$latlon_g----
# -------------------------------------------------------------------------

knots_latlong <- as.data.frame(fit_f$spatial_list$latlon_g)
knots <- c(1:100)
knots <- as.data.frame(knots)
knots_latlong <- cbind(knots_latlong,knots)
knots_sf <- st_as_sf(knots_latlong, coords = c('Lon', 'Lat'), crs = "+proj=longlat +datum=WGS84")
concav_survey2 <-st_as_sf(survey, coords = c("Lon_i", "Lat_i"),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>% 
  concaveman::concaveman(concavity = 2)

w <- ggplot() +
 # geom_sf(data = concav_survey2,fill=NA, size=2,color="pink")+
  geom_sf(data = knots_sf,color="red")+
  geom_sf(data = shapefile,fill=NA,color="blue") + 
  ggtitle("EBS") +
  geom_point(survey,mapping=aes(x=Lon_i,y=Lat_i),color="blue")+
  geom_sf(data = world,fill="black",color=NA) + coord_sf(xlim=c(-180,-155), ylim=c(54,63.5),expand = FALSE)



png(paste('EBSarea','.png',sep=''), height =14 , width = 20, units = 'cm', res=500)#, res=500)
w
dev.off()


# Check data vs knot distribution
w+ geom_point(fit_f$data_frame,mapping=aes(x=Lon_i,y=Lat_i),color="green3",size=2)
u <- ggplot() +
  geom_sf(data = knots_sf,color="red")+
  geom_point(fit_f$data_frame,mapping=aes(x=Lon_i,y=Lat_i),color="green3",size=2) 
u

# -------------------------------------------------------------------------
# 3-  Make an indicator Inside_g representing whether each grid cell is within (Inside_g[g] = 1) ---------
# the EBS survey or not (Inside_g[g] = 0)
# -------------------------------------------------------------------------

dim(knots_sf)
pnts <- st_difference( knots_sf,shapefile)
pnts <- st_transform(pnts, "+proj=longlat +datum=WGS84")
q <- ggplot() +
#  geom_sf(data = concav_survey2,fill=NA, size=2,color="pink")+
  geom_sf(data = pnts,color="red")+
  geom_sf(data = shapefile,fill=NA,color="blue") + 
  ggtitle("EBS") +
  geom_point(survey,mapping=aes(x=Lon_i,y=Lat_i),color="blue")+
  geom_sf(data = world,fill="black",color=NA) + coord_sf(xlim=c(-180,-155), ylim=c(54,63.5),expand = FALSE)
q


pnts_coord <- cbind(st_coordinates(pnts), pnts$knots)
index <- as.data.frame(rep(0,dim(pnts_coord)[1]))
colnames(index)=c('index')
pnts_coord <- cbind(pnts_coord, index)
colnames(pnts_coord) <- c("Lon", "Lat","knots","index")

Index_g <- left_join(knots_latlong,pnts_coord)
Index_g <- Index_g %>% replace(is.na(.),1)

# - 4 Get the estimate of biomass fit$Report$d_gct, and average across seasons u and conditions v to get d_guv --------
# -------------------------------------------------------------------------

# -- fisheries CPUE --------------------------------------------------------
load("04_outputs/user_grid/EnvCxSeason/c(3,3,3,3,2,2)_grid1/fit.RData")
fit_f <- fit
dim(fit_f$Report$D_gct)
fit$data_list$n_i
fit_f$Report$Phi1_ik

fit_f$data_frame

#plot data points
q + geom_point(fit_f$data_frame,mapping=aes(x=Lon_i,y=Lat_i),color="green3")


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
      D_cold_gut[g,u,t] <- fit_f$Report$D_gct[g,1,cold[t]] * 
        # a_k[g,1] *
        exp(fit_f$Report$Phi1_sk[g,cold_season[u]]) 
    }
  }
}

for ( t in 1:length(warm)){
  for ( u in 1 :length(warm_season)){
    for ( g in 1:n_knot){
      D_warm_gut[g,u,t] <- fit_f$Report$D_gct[g,1,warm[t]] *
        #a_k[g,1] * 
        exp(fit_f$Report$Phi1_sk[g,warm_season[u]])  
    }
  }
}
dim(fit_f$Report$D_gct)

# - Average accross all years ---------------------------------------------
D_cold_gu <- apply(D_cold_gut,c(1,2),mean)
D_warm_gu <- apply(D_warm_gut,c(1,2),mean)


# I_uv: Total abundance acrross space and season and env conditions -------------
I_coldu <- apply(D_cold_gu,2,sum)
I_warmu <- apply(D_warm_gu,2,sum)
I_uv <- rbind(I_warmu,I_coldu)
colnames(I_uv)<-  c("Early","Int","Late")
rownames(I_uv) <- c("Warm","Cold")

# P_uv : abundance within survey ------------------------------------------
# -- cold
colnames(D_cold_gu) <- c("Early","Int","Late")
D_cold_gu_temp <- cbind(D_cold_gu,Index_g)

P_cold_gu <- D_cold_gu_temp %>% mutate(Early= Early*index, Int=Int*index, Late=Late*index)
P_cold_u <- P_cold_gu %>% mutate(Early=sum(Early), Int=sum(Int), Late=sum(Late)) 
P_cold_u <- P_cold_u[1,c(1:3)]

# -- warm
colnames(D_warm_gu) <- c("Early","Int","Late")
D_warm_gu_temp <- cbind(D_warm_gu,Index_g)

P_warm_gu <- D_warm_gu_temp %>% mutate(Early= Early*index, Int=Int*index, Late=Late*index)
P_warm_u <- P_warm_gu %>% mutate(Early=sum(Early), Int=sum(Int), Late=sum(Late)) 
P_warm_u <- P_warm_u[1,c(1:3)]

# P_u : merge
P_uv <- rbind(P_warm_u,P_cold_u)


# Calcul O_uv -------------------------------------------------------------
O_uv <- P_uv/I_uv
rownames(O_uv) <- c("warm","cold")



# plot --------------------------------------------------------------------

O_uv <- melt(O_uv)
O_uv <- cbind(O_uv, c("warm","cold","warm","cold","warm","cold"))
O_uv <- O_uv %>% mutate(model="Density_tot")
colnames(O_uv) <- c("Season", "Overlap", "EnvCond","model")

ggplot() + geom_point(O_uv, mapping=aes(x=Season, y=Overlap, color=EnvCond),size=4)+
  scale_color_manual(values=c( "#1E90FF","#FF0000"))



# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# CALCULATE ONLY WITH OMEGA AND PHI ---------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------


dim(fit_f$Report$Phi1_sk)
D_cold_gu <- array(NA,dim=c(n_knot,length(cold_season)))
D_warm_gu <- array(NA,dim=c(n_knot,length(warm_season)))


for ( u in 1 :length(cold_season)){
  for ( g in 1:n_knot){
    D_cold_gu[g,u] <- exp(fit_f$Report$Omega1_gc[g,1])*
      #a_k[g,1] * 
      exp(fit_f$Report$Phi1_sk[g,cold_season[u]]) 
  }
}


for ( u in 1 :length(warm_season)){
  for ( g in 1:n_knot){
    D_warm_gu[g,u] <- exp(fit_f$Report$Omega1_gc[g,1])*
      #a_k[g,1]* 
      exp(fit_f$Report$Phi1_sk[g,warm_season[u]])  
  }
}


# I_uv: Total abundance acrross space and season and env conditions -------------
I_coldu <- apply(D_cold_gu,2,sum)
I_warmu <- apply(D_warm_gu,2,sum)
I_uv <- rbind(I_warmu,I_coldu)
colnames(I_uv)<-  c("Early","Int","Late")
rownames(I_uv) <- c("Warm","Cold")

# P_uv : abundance within survey ------------------------------------------
# -- cold
colnames(D_cold_gu) <- c("Early","Int","Late")
D_cold_gu_temp <- cbind(D_cold_gu,Index_g)

P_cold_gu <- D_cold_gu_temp %>% mutate(Early= Early*index, Int=Int*index, Late=Late*index)
P_cold_u <- P_cold_gu %>% mutate(Early=sum(Early), Int=sum(Int), Late=sum(Late)) 
P_cold_u <- P_cold_u[1,c(1:3)]

# -- warm
colnames(D_warm_gu) <- c("Early","Int","Late")
D_warm_gu_temp <- cbind(D_warm_gu,Index_g)

P_warm_gu <- D_warm_gu_temp %>% mutate(Early= Early*index, Int=Int*index, Late=Late*index)
P_warm_u <- P_warm_gu %>% mutate(Early=sum(Early), Int=sum(Int), Late=sum(Late)) 
P_warm_u <- P_warm_u[1,c(1:3)]

# P_u : merge
P_uv <- rbind(P_warm_u,P_cold_u)


# Calcul O_uv -------------------------------------------------------------
O_uv_omega <- P_uv/I_uv
rownames(O_uv_omega ) <- c("warm","cold")

O_uv_omega <- melt(O_uv_omega)
O_uv_omega <- cbind(O_uv_omega, c("warm","cold","warm","cold","warm","cold"))
colnames(O_uv_omega) <- c("Season", "Overlap", "EnvCond")
#O_uv_omega$Season <- as.numeric(O_uv_omega$Season)
O_uv_omega <- O_uv_omega %>% mutate()
O_uv_omega <- O_uv_omega %>% mutate(model="Density_expected")

Overlap <- rbind(O_uv_omega,O_uv)


p<-ggplot(Overlap, mapping=aes(x=Season, y=Overlap, color=EnvCond,shape=model)) +
  scale_color_manual(values=c( "#1E90FF","#FF0000"))+   geom_point(size=3)#+geom_line()
getwd()
png(paste('Overlap_comp','.png',sep=''), height =14 , width = 20, units = 'cm', res=500)#, res=500)
p
dev.off()


# -------------------------------------------------------------------------
# densiites
# ------------------------------------------------------------------------
# -------------------------------------------------------------------------


dim(fit_f$Report$Phi1_sk)
D_cold_gu <- array(NA,dim=c(n_knot,length(cold_season)))
D_warm_gu <- array(NA,dim=c(n_knot,length(warm_season)))


for ( u in 1 :length(cold_season)){
  for ( g in 1:n_knot){
    D_cold_gu[g,u] <- exp(fit_f$Report$Omega1_gc[g,1]) * exp(fit_f$Report$Phi1_sk[g,cold_season[u]] + 
                                                               fit_f$ParHat$lambda1_k[cold_season[u]] ) 
  }
}

#fit_f$ParHat$lambda1_k


for ( u in 1 :length(warm_season)){
  for ( g in 1:n_knot){
    D_warm_gu[g,u] <- exp(fit_f$Report$Omega1_gc[g,1])* exp(fit_f$Report$Phi1_sk[g,warm_season[u]]+
                                                              fit_f$ParHat$lambda1_k[warm_season[u]])  
  }
}


# I_uv: Total abundance acrross space and season and env conditions -------------
I_coldu <- apply(D_cold_gu,2,sum)
I_warmu <- apply(D_warm_gu,2,sum)
I_uv <- rbind(I_warmu,I_coldu)
colnames(I_uv)<-  c("Early","Int","Late")
rownames(I_uv) <- c("Warm","Cold")

# P_uv : abundance within survey ------------------------------------------
# -- cold
colnames(D_cold_gu) <- c("Early","Int","Late")
D_cold_gu_temp <- cbind(D_cold_gu,Index_g)

P_cold_gu <- D_cold_gu_temp %>% mutate(Early= Early*index, Int=Int*index, Late=Late*index)
P_cold_u <- P_cold_gu %>% mutate(Early=sum(Early), Int=sum(Int), Late=sum(Late)) 
P_cold_u <- P_cold_u[1,c(1:3)]

# -- warm
colnames(D_warm_gu) <- c("Early","Int","Late")
D_warm_gu_temp <- cbind(D_warm_gu,Index_g)

P_warm_gu <- D_warm_gu_temp %>% mutate(Early= Early*index, Int=Int*index, Late=Late*index)
P_warm_u <- P_warm_gu %>% mutate(Early=sum(Early), Int=sum(Int), Late=sum(Late)) 
P_warm_u <- P_warm_u[1,c(1:3)]

# P_u : merge
P_uv <- rbind(P_warm_u,P_cold_u)


# Calcul O_uv -------------------------------------------------------------
O_uv_omega_d <- P_uv/I_uv
rownames(O_uv_omega_d ) <- c("warm","cold")

