rm(list=ls())

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
library(patchwork)
# Fisheries ---------------------------------------------------------------
work_dir <- "C:/Users/Maxime/Documents/Git/Flat_fish_2021/04_outputs/user_grid_outEBS/no_sex/finescale/c(3,3,3,3,2,2)/no2000/av_season/"
#"C:/Users/Maxime/Documents/Git/Flat_fish_2021/04_outputs/user_grid_outEBS/sex/female/no2000/"
load(paste0(work_dir,"fit.RData"))
fit_f <- fit
dim(fit_f$Report$D_gct)
fit_f$data_frame
getwd()
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Spatial features --------------------------------------------------------
# -------------------------------------------------------------------------

# extra grid
Region = "Eastern_Bering_Sea"
MapDetails_List <- FishStatsUtils:: make_map_info( "Region"=Region,spatial_list = fit$spatial_list, "NN_Extrap"= fit$spatial_lis$PolygonList$NN_Extrap,                                                     
                                                   "Extrapolation_List"=fit$extrapolation_list ,fine_scale = TRUE)  
PlotDF1 = subset(MapDetails_List$PlotDF)
Knot_to_grid <- PlotDF1[,c(1,2,3)]
colnames(Knot_to_grid) <- c("Lat","Lon","cell")

# Shape file
shape_path <- "C:/Users/Maxime/Documents/POST_DOC_SNOWCRAB/Spatial_snow_crab/outputs/function_r/simulator/"
coast_shapefile <- paste(shape_path, "ne_50m_land.shp", sep="")
ocean <- readOGR(coast_shapefile)

# survey shapefile
getwd()
shapefile <- st_read("C:/Users/Maxime/Documents/Git/Spatial_snow_crab_2021/02_transformed_data/EBSshelf/EBSshelf.shp")
crs_ref <- "+proj=longlat +ellps=WGS84  +datum=WGS84 +no_defs"
EBS <- st_transform(shapefile,crs=crs_ref)

# world
world <- ne_countries(scale = "medium", returnclass = "sf")

# -------------------------------------------------------------------------
# Generates tibble output -------------------------------------------------
# -------------------------------------------------------------------------

# Demsity------------------------------------------------------------------
# -------------------------------------------------------------------------
# warm
dim(D_warm_gut)
D_gut_warm <- array(NA, dim=dim(D_warm_gut[,,]),dimnames=list(c(1:fit_f$spatial_list$n_g),c("Early","Int.","Late"),warm_year)) 
for (k in 1:fit_f$spatial_list$n_g){
  D_gut_warm[k,,] <- D_warm_gut[k,,]
}

D_gut_warm <- melt(D_gut_warm)
colnames(D_gut_warm) <- c("cell","Season","Years","value")

D_gut_warm_grid <- left_join(as.data.frame(Knot_to_grid),as.data.frame(D_gut_warm))%>% mutate(En.Cond = "Warm")


#cold
D_gut_cold <- array(NA, dim=dim(D_cold_gut[,,]),dimnames=list(c(1:fit_f$spatial_list$n_g),c("Early","Int.","Late"),cold_year)) 
for (k in 1:fit_f$spatial_list$n_g){
  D_gut_cold[k,,] <- D_cold_gut[k,,]
}

D_gut_cold <- melt(D_gut_cold)
colnames(D_gut_cold) <- c("cell","Season","Years","value")
D_gut_cold_grid <- left_join(as.data.frame(Knot_to_grid),as.data.frame(D_gut_cold)) %>% mutate(En.Cond = "Cold")

D_gut_grid <- bind_rows(D_gut_cold_grid,D_gut_warm_grid)

season

xlims <- c(-179.5,-155) #range(pretty(D_gut_grid$Lon))
ylims <- range(pretty(D_gut_grid$Lat))

for ( i in 1 : length(season)){
  
 map1 <- ggplot()+ 
    geom_point(data=D_gut_grid %>% filter(Season== season[i] ),aes(x=Lon,y=Lat,col=log(value)),size=0.5)+
    #geom_sf(data=Ab_map %>% filter(size_class==1),aes(col=log(A_pkt)))+
    scale_color_viridis() + labs(colour= "Predicted biomass")+
 #   geom_sf(data=EBS,fill=NA, color="brown")+
    geom_sf(data = sebs_layers$bathymetry,size=0.2) +
    geom_sf(data=world, col=NA, fill="black")+
    coord_sf(xlim = xlims, ylim = ylims)+
    facet_wrap( ~ Years)+
    theme_bw(base_size = 8)

 map1 <- ggplot()+ 
   geom_point(data=D_gut_grid %>% filter(Season== season[i] ),aes(x=Lon,y=Lat,col=log(value)),size=1)+
   #geom_sf(data=Ab_map %>% filter(size_class==1),aes(col=log(A_pkt)))+
   scale_color_viridis() + labs(colour= "Predicted biomass")+
   #   geom_sf(data=EBS,fill=NA, color="brown")+
   geom_sf(data = sebs_layers$bathymetry,size=0.2) +
   geom_sf(data=world, col=NA, fill="black")+
   coord_sf(xlim = xlims, ylim = ylims)+
   facet_wrap( ~ Years)
}  # theme_bw(base_size = 8)
 
 
 map1 <- ggplot()+ 
   geom_point(data=D_gut_grid %>% filter(Years %in% c(2001:2010)),aes(x=Lon,y=Lat,col=log(value)),size=1)+
   #geom_sf(data=Ab_map %>% filter(size_class==1),aes(col=log(A_pkt)))+
   scale_color_viridis() + labs(colour= "Predicted biomass")+
   #   geom_sf(data=EBS,fill=NA, color="brown")+
   geom_sf(data = sebs_layers$bathymetry,size=0.2) +
   geom_sf(data=world, col=NA, fill="black")+
   coord_sf(xlim = xlims, ylim = ylims)+
   facet_grid(Years ~ Season)+
   theme(axis.text = element_blank())
  # theme_bw(base_size = 8)
 
  
  ggsave((paste0(work_dir, "map1.png")),plot=map1,
         width = 27,
         height = 30,
         units = "cm")
  
  map2 <- ggplot()+ 
    geom_point(data=D_gut_grid %>% filter(Years %in% c(2011:2019)),aes(x=Lon,y=Lat,col=log(value)),size=1)+
    #geom_sf(data=Ab_map %>% filter(size_class==1),aes(col=log(A_pkt)))+
    scale_color_viridis() + labs(colour= "Predicted biomass")+
    #   geom_sf(data=EBS,fill=NA, color="brown")+
    geom_sf(data = sebs_layers$bathymetry,size=0.2) +
    geom_sf(data=world, col=NA, fill="black")+
    coord_sf(xlim = xlims, ylim = ylims)+
    facet_grid(Years ~ Season)+
    theme(axis.text = element_blank())
    #theme_bw(base_size = 8)
  
  
  ggsave((paste0(work_dir, "map2.png")),plot=map2,
         width = 27,
         height = 30,
         units = "cm")



# Spatial effect-----------------------------------------------------------
# -------------------------------------------------------------------------
# warm

phi <- (fit_f$Report$Phi1_gk)
phi1 <- as.data.frame(cbind(phi[,1],"Cold","Early",Knot_to_grid ))
phi2 <- as.data.frame(cbind(phi[,2],"Warm","Early",Knot_to_grid))
phi3 <- as.data.frame(cbind(phi[,3],"Cold","Int.",Knot_to_grid))
phi4 <- as.data.frame(cbind(phi[,4],"Warm","Int.",Knot_to_grid))
phi5 <- as.data.frame(cbind(phi[,5],"Cold","Late",Knot_to_grid))
phi6 <- as.data.frame(cbind(phi[,6],"Warm","Late",Knot_to_grid))

phi_tibble <- rbind(phi1,phi2,phi3,phi4,phi5,phi6)

colnames(phi_tibble) <- c("value","EnvCond","Seasons","Lat","Lon","cell")
phi_tibble$value <- as.numeric(phi_tibble$value)
phi_tibble$Lat <- as.numeric(phi_tibble$Lat)
phi_tibble$Lon <- as.numeric(phi_tibble$Lon)

xlims <- range(pretty(phi_tibble$Lon))
ylims <- range(pretty(phi_tibble$Lat))


map1 <- ggplot()+ 
  geom_point(data=as_tibble(phi_tibble),aes(x=Lon,y=Lat,col=value),size=1.2)+
  #geom_sf(data=Ab_map %>% filter(size_class==1),aes(col=log(A_pkt)))+
  scale_color_viridis() + labs(colour= "Spatial effect")+
#  geom_sf(data=EBS,fill=NA, color="brown")+
  geom_sf(data = sebs_layers$bathymetry,size=0.2) +
  geom_sf(data=world, col=NA, fill="black")+
  coord_sf(xlim = xlims, ylim = ylims)+
  facet_grid( EnvCond ~ Seasons)
  #theme_bw(base_size = 8)

ggsave((paste0(work_dir, 'Spatialeffect.png')),
       width = 25,
       height = 20,
       units = "cm")


# -------------------------------------------------------------------------
#Significant effect
# -------------------------------------------------------------------------


load(paste0(work_dir,"/fit_sign.RData"))


# Re sample  --------------------------------------------------------------
sample_true = sample_variable( Sdreport=fit_sign$parameter_estimates$SD,
                               Obj=fit_sign$tmb_list$Obj, variable_name="Phi1_gk" )
sample_false = sample_variable( Sdreport=fit_sign$parameter_estimates$SD,n_samples =250,
                                Obj=fit_sign$tmb_list$Obj, variable_name="Phi1_gk",sample_fixed=FALSE )

# TRUE
var <- apply(sample_true,c(1,2),var)
# FALSE-> work better with false
var <- apply(sample_false,c(1,2),var)
sd <- sqrt(var)
Wt <- fit$Report$Phi1_gk/sd # t value from t test in summary

# calcul pvalue : same as pvalue in summary() test
pvalue = pchisq(Wt^2, df = 1, lower = F) #
pvalue2 =pnorm(-abs(Wt), 0, 1, lower = T) + pnorm(abs(Wt), 0, 1, lower = F) #
round(pvalue-pvalue2, 5)

pvaluenNA <- ifelse(pvalue>0.05, pvalue==NA, pvalue)


# Map  --------------------------------------------------------------------

pvalue <- pvaluenNA
pvalue1 <- as.data.frame(cbind(pvalue[,1],"Cold","Early",Knot_to_grid ))
pvalue2 <- as.data.frame(cbind(pvalue[,2],"Warm","Early",Knot_to_grid))
pvalue3 <- as.data.frame(cbind(pvalue[,3],"Cold","Int.",Knot_to_grid))
pvalue4 <- as.data.frame(cbind(pvalue[,4],"Warm","Int.",Knot_to_grid))
pvalue5 <- as.data.frame(cbind(pvalue[,5],"Cold","Late",Knot_to_grid))
pvalue6 <- as.data.frame(cbind(pvalue[,6],"Warm","Late",Knot_to_grid))

pavlue_tibble <- rbind(pvalue1,pvalue2,pvalue3,pvalue4,pvalue5,pvalue6)

colnames(pavlue_tibble) <- c("value","EnvCond","Seasons","Lat","Lon","cell")
pavlue_tibble$value <- as.numeric(pavlue_tibble$value)
pavlue_tibble$Lat <- as.numeric(pavlue_tibble$Lat)
pavlue_tibble$Lon <- as.numeric(pavlue_tibble$Lon)

xlims <- range(pretty(pavlue_tibble$Lon))
ylims <- range(pretty(pavlue_tibble$Lat))


## bathy
library(akgfmaps)
sebs_layers <- akgfmaps::get_base_layers(select.region = "sebs",
                                         set.crs =crs_ref)

map1 <- ggplot()+ 
  geom_point(data=as_tibble(pavlue_tibble),aes(x=Lon,y=Lat,col=value),size=1.2)+
  #geom_sf(data=Ab_map %>% filter(size_class==1),aes(col=log(A_pkt)))+
 # geom_sf(data=EBS,fill=NA, color="brown")+
  geom_sf(data = sebs_layers$bathymetry,size=0.2) + 
  scale_color_viridis(direction=-1,  na.value = "grey94") + labs(colour= "Spatial effect")+
  geom_sf(data=world, col=NA, fill="black")+
  coord_sf(xlim = xlims, ylim = ylims)+
  facet_grid( EnvCond ~ Seasons)+
  theme_bw(base_size = 8)

ggsave((paste0(work_dir, 'Pavlue_sign.png')),
       width = 25,
       height = 20,
       units = "cm")


# -------------------------------------------------------------------------
phi_tibble2 = phi_tibble

for (i in 1:dim(phi_tibble)[1]){
if (is.na(pavlue_tibble$value[i])==TRUE) {phi_tibble2$value[i]=NA}
}


map1 <- ggplot()+ 
  geom_point(data=as_tibble(phi_tibble2),aes(x=Lon,y=Lat,col=value),size=1.2)+
  #geom_sf(data=Ab_map %>% filter(size_class==1),aes(col=log(A_pkt)))+
  # geom_sf(data=EBS,fill=NA, color="brown")+
  geom_sf(data = sebs_layers$bathymetry,size=0.2) + 
  scale_color_viridis(  na.value = "grey94") + labs(colour= "Spatial effect")+
  geom_sf(data=world, col=NA, fill="black")+
  coord_sf(xlim = xlims, ylim = ylims)+
  facet_grid( EnvCond ~ Seasons)+
  theme_bw(base_size = 8)

ggsave((paste0(work_dir, 'Spatialeffect_sign.png')),
       width = 25,
       height = 20,
       units = "cm")

