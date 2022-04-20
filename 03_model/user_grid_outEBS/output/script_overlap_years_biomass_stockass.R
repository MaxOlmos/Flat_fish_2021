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
work_dir <- "C:/Users/Maxime/Documents/Git/Flat_fish_2021/04_outputs/user_grid_outEBS/no_sex/finescale/c(3,3,3,3,2,2)/no2000/sum/av_season/"
#"C:/Users/Maxime/Documents/Git/Flat_fish_2021/04_outputs/user_grid_outEBS/sex/male/no2000/"
load(paste0(work_dir,"fit.RData"))
fit_f <- fit
dim(fit_f$Report$D_gct)
fit_f$data_frame
getwd()
# Overlap steps -----------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# - 1 Get the shapefile for the EBS survey e.g. here:  --------------------
#https://github.com/James-Thorson-NOAA/FishStatsUtils/tree/main/inst/region_shapefiles/EBSshelf
# -------------------------------------------------------------------------
# world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# survey shapefile
shapefile <- st_read("02_transformed_data/EBSshelf/EBSshelf.shp")
crs_ref <- "+proj=longlat +ellps=WGS84  +datum=WGS84 +no_defs"
EBS <- st_transform(shapefile,crs=crs_ref)


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

sum(as.data.frame(fit_f$data_list$a_gl) - as.numeric(fit_f$extrapolation_list$a_el))


extgrid_latlong <- cbind(as.data.frame(fit_f$spatial_list$latlon_g),as.data.frame(fit_f$data_list$a_gl))
colnames(extgrid_latlong) <- c("Lat","Lon","area_grid")
head(extgrid_latlong)
extgrid <- c(1:dim(extgrid_latlong)[1])
extgrid <- as.data.frame(extgrid)
extgrid_latlong <- cbind(extgrid_latlong,extgrid)
extgrid_sf <- st_as_sf(extgrid_latlong, coords = c('Lon', 'Lat'), crs = crs_ref)

w <- ggplot() +
  # geom_sf(data = concav_survey2,fill=NA, size=2,color="pink")+
  geom_sf(data = extgrid_sf,color="red")+
  geom_sf(data = shapefile,fill=NA,color="blue") + 
  ggtitle("EBS") +
  geom_point(survey,mapping=aes(x=Lon_i,y=Lat_i),color="blue")+
  geom_sf(data = world,fill="black",color=NA) + coord_sf(xlim=c(-180,-155), ylim=c(54,63.5),expand = FALSE)



png(paste0(work_dir,'EBSarea','.png'), height =14 , width = 20, units = 'cm', res=500)#, res=500)
w
dev.off()


# Check data vs knot distribution
w+ geom_point(fit_f$data_frame,mapping=aes(x=Lon_i,y=Lat_i),color="green3",size=2)
u <- ggplot() +
  geom_sf(data = extgrid_sf,color="red")+
  geom_point(fit_f$data_frame,mapping=aes(x=Lon_i,y=Lat_i),color="green3",size=2) 
u

# -------------------------------------------------------------------------
# 3-  Make an indicator Inside_g representing whether each grid cell is within (Inside_g[g] = 1) ---------
# the EBS survey or not (Inside_g[g] = 0)
# -------------------------------------------------------------------------

dim(extgrid_sf)
pnts <- st_difference( extgrid_sf,EBS)
pnts <- st_transform(pnts, "+proj=longlat +datum=WGS84")
q <- ggplot() +
  #  geom_sf(data = concav_survey2,fill=NA, size=2,color="pink")+
  geom_sf(data = pnts,color="red")+
  geom_sf(data = shapefile,fill=NA,color="blue") + 
  ggtitle("EBS") +
  geom_point(survey,mapping=aes(x=Lon_i,y=Lat_i),color="blue")+
  geom_sf(data = world,fill="black",color=NA) + coord_sf(xlim=c(-180,-155), ylim=c(54,63.5),expand = FALSE)
q

pnts_coord <- cbind(st_coordinates(pnts), pnts$extgrid)
index <- as.data.frame(rep(0,dim(pnts_coord)[1]))
colnames(index)=c('index')
pnts_coord <- cbind(pnts_coord, index)
colnames(pnts_coord) <- c("Lon", "Lat","extgrid","index")

Index_g <- left_join(extgrid_latlong,as_tibble(pnts_coord) %>% dplyr::select("extgrid","index"))
Index_g <- Index_g %>% replace(is.na(.),1)

ggplot() +
  geom_point(data=Index_g, aes(x=Lon,y=Lat,color=index))+
  # geom_point(data=pnts_coord, aes(x=Lon,y=Lat),color="red",shape=3)+
  geom_sf(data = EBS,fill=NA,color="red")#+
#  geom_sf(data=extgrid_sf,fill=NA, color="green")

# - 4 Get the estimate of biomass fit$Report$d_gct, and average across seasons u and conditions v to get d_guv --------
# -------------------------------------------------------------------------

# -- fisheries CPUE --------------------------------------------------------
#plot data points
q + geom_point(fit_f$data_frame,mapping=aes(x=Lon_i,y=Lat_i),color="green3")


# Calculate abundance because we want proportion
# -------------------------------------------------------------------------
fit_f$data_frame
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
D_cold_gut <- array(NA,dim=c(dim(fit_f$Report$Phi1_gk)[1],length(cold_season),length(cold)))
D_warm_gut <- array(NA,dim=c(dim(fit_f$Report$Phi1_gk)[1],length(warm_season),length(warm)))

for ( t in 1:length(cold)){
  for ( u in 1 :length(cold_season)){
    for ( g in 1:dim(fit_f$Report$Phi1_gk)[1]){
      D_cold_gut[g,u,t] <- fit_f$Report$D_gct[g,,cold[t]] * Index_g$area_grid[g] *exp(fit_f$Report$Phi1_gk[g,cold_season[u]]) 
    }
  }
}

for ( t in 1:length(warm)){
  for ( u in 1 :length(warm_season)){
    for ( g in 1:dim(fit_f$Report$Phi1_gk)[1]){
      D_warm_gut[g,u,t] <- fit_f$Report$D_gct[g,1,warm[t]] * Index_g$area_grid[g] * exp(fit_f$Report$Phi1_gk[g,warm_season[u]])  
    }
  }
}

# -------------------------------------------------------------------------
# Calculate overlap per year-----------------------------------------------
# -------------------------------------------------------------------------

dim(D_warm_gut)
dim(D_cold_gut)

# I_uv: Total abundance acrross  season and env conditions -------------

#cold
I_cold_ut <- apply(D_cold_gut,c(2,3),sum)
colnames(I_cold_ut) <- cold_year
rownames(I_cold_ut) <- season
I_cold_ut <- melt(I_cold_ut)
colnames(I_cold_ut) <- c("Seasons","Years","Ivalue")
I_cold_ut <- I_cold_ut %>% mutate(Env.Cond = "Cold")

#warm
I_warm_ut <- apply(D_warm_gut,c(2,3),sum)
colnames(I_warm_ut) <- warm_year
rownames(I_warm_ut) <- season
I_warm_ut <- melt(I_warm_ut)
colnames(I_warm_ut) <- c("Seasons","Years","Ivalue")
I_warm_ut <- I_warm_ut %>% mutate(Env.Cond = "Warm")


I_uvt <- rbind(I_warm_ut,I_cold_ut)


# P_uv : abundance within survey ------------------------------------------

# -- cold
dim(D_cold_gut)
D_cold_gut <- array(D_cold_gut,dim=dim(D_cold_gut), dimnames=list(c(1:dim(D_cold_gut)[1]),season,cold_year))
D_cold_gut_tibble <- melt(D_cold_gut)
colnames(D_cold_gut_tibble) <- c("extgrid","Seasons","Years","value")


D_cold_gut_temp <- left_join(D_cold_gut_tibble,Index_g)
P_cold_gut_index <- D_cold_gut_temp %>% mutate(value_P= value*index) %>%
  dplyr::group_by(Years,Seasons) %>% 
  dplyr::summarise(Pvalue= sum(value_P)) %>% 
  dplyr::mutate(Env.Cond = "Cold")


# warm  
D_warm_gut <- array(D_warm_gut,dim=dim(D_warm_gut), dimnames=list(c(1:dim(D_warm_gut)[1]),season,warm_year))
D_warm_gut_tibble <- melt(D_warm_gut)
colnames(D_warm_gut_tibble) <- c("extgrid","Seasons","Years","value")


D_warm_gut_temp <- left_join(D_warm_gut_tibble,Index_g)
P_warm_gut_index <- D_warm_gut_temp %>% mutate(value_P= value*index) %>%
  dplyr::group_by(Years,Seasons) %>% 
  dplyr::summarise(Pvalue= sum(value_P)) %>% 
  mutate(Env.Cond = "Warm")


# P_u : merge
P_uvt <- rbind(P_warm_gut_index,P_cold_gut_index)

# Calcul O_uvt -------------------------------------------------------------
O_uvt <- left_join(P_uvt,as_tibble(I_uvt))

O_uvt <- O_uvt %>% mutate(Overlap=Pvalue/Ivalue)

# calcul O_uv -------------------------------------------------------------

O_uv <- O_uvt %>% dplyr::group_by(Seasons,Env.Cond) %>% dplyr::summarise(Overlap_m=mean(Overlap))


# plot --------------------------------------------------------------------
O_uvt_2 <- left_join(O_uvt,O_uv)

# calculate sd
O_uvt_3 <- O_uvt_2 %>% dplyr::group_by(Seasons, Env.Cond) %>% dplyr::summarise(var=var(Overlap),sd=100*sd(Overlap))
O_uvt_4 <- left_join(O_uvt_2,O_uvt_3)

p <- ggplot() + 
  geom_line(O_uvt_4, mapping=aes(x=Years, y=Overlap_m,color=Env.Cond,size=var),alpha=0.2)+
  ylab("Overlap")+
  geom_line(O_uvt_4, mapping=aes(x=Years, y=Overlap),size=1)+
  geom_point(O_uvt, mapping=aes(x=Years, y=Overlap, color=Env.Cond),size=2)+
  scale_color_manual(values=c( "#1E90FF","#FF0000"))+ guides(size = FALSE)

pseason <- p + facet_wrap( ~ Seasons) + ggtitle('a)')
pseason


#O_uvt_male <- O_uvt_4 %>% mutate(sex="male")
save(O_uvt,file=paste0(work_dir, "O_uvt.RData"))


# Catchabiltity from stock assessment
# -------------------------------------------------------------------------

SA_outputs <- read.table("C:/Users/Maxime/Documents/Git/Flat_fish_2021/01_data/stock_assessment_outputs/outputs.csv", header=T, sep=",")
SA_outputs <- as_tibble(SA_outputs)
colnames(SA_outputs)[1] <- c("Years")
colnames(SA_outputs)[2] <- c("q")

SA_outputs <- SA_outputs %>% dplyr::filter(Years %in% c(2001:2019))

Catchability <- SA_outputs %>% dplyr::select(Years, q) %>% mutate(Type="Catchability")
Catchability <- Catchability %>% 
  mutate(Index = (q-mean(q))/sd(q)) %>% dplyr::select(Years,Index,Type)
# ------------------------------------------------

O_index_int <- O_uvt %>% filter(Seasons =="Int.")
O_index <- (O_index_int$Overlap - mean(O_index_int$Overlap))/sd(O_index_int$Overlap) 
O_index <- as.data.frame(O_index)

O_index_int <- cbind(O_index,O_index_int$Years)
colnames(O_index_int) <- c("Index", "Years")
O_index_int <- O_index_int %>% mutate(Type = "Overlap")

Index_comp <- rbind(O_index_int,Catchability %>% mutate(Years=Years+1))
Index_comp <- rbind(O_index_int,Catchability)

Overlap <-arrange(O_index_int,Years) %>% dplyr::select(Years, Index,Type)
getwd()
save(Overlap,file=paste0(work_dir, "Overlap.RData"))

# no smooth
#----------------
pcompare <- ggplot(Index_comp, mapping=aes(x=Years, y=Index,color=Type),size=2) + geom_point(size=2)+ geom_line(size=1.5)+
  #geom_smooth(aes(color=Type,fill=Type))+
  #geom_point(Index_comp, mapping=aes(x=Years, y=IA_index, color=Env.Cond),size=4)+
  scale_color_manual(values=c( "yellow3","purple"))
p_comp_enc <- pcompare +  
  geom_rect(aes(xmin = 2001, xmax = 2006, ymin = -Inf, ymax = Inf),
            fill = "pink",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2006, xmax = 2014, ymin = -Inf, ymax = Inf),
            fill = "lightblue",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2014, xmax = 2017, ymin = -Inf, ymax = Inf),
            fill = "pink",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2017, xmax = 2018, ymin = -Inf, ymax = Inf),
            fill = "lightblue",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2018, xmax = 2019, ymin = -Inf, ymax = Inf),
            fill = "pink",col=NA, alpha = 0.01) + ggtitle('b)')


ggsave(paste0(work_dir,'Catchability_vs_overalap.png'), p_comp_enc)

# smooth
#--------------
pcompare_s <- ggplot(Index_comp, mapping=aes(x=Years, y=Index,color=Type),size=2) + geom_point(size=2)+
  geom_smooth(aes(color=Type),se=FALSE,size=2)+
  scale_color_manual(values=c( "purple","yellow3"))+
  scale_fill_manual(values = c("purple", "yellow3"))  
pcompare_s


pcompare_s_enc <- pcompare_s +  
  geom_rect(aes(xmin = 2001, xmax = 2006, ymin = -Inf, ymax = Inf),
            fill = "pink",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2006, xmax = 2014, ymin = -Inf, ymax = Inf),
            fill = "lightblue",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2014, xmax = 2017, ymin = -Inf, ymax = Inf),
            fill = "pink",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2017, xmax = 2018, ymin = -Inf, ymax = Inf),
            fill = "lightblue",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2018, xmax = 2019, ymin = -Inf, ymax = Inf),
            fill = "pink",col=NA, alpha = 0.01) + ggtitle('c)')



pseason


p <- pseason  / (p_comp_enc + pcompare_s_enc)
ggsave(paste0(work_dir,'Overlap_ab_catchability.png'), p)
getwd()



# -------------------------------------------------------------------------
# calculate residual  -----------------------------------------------------
# -------------------------------------------------------------------------

# residuals from stock assessment
# -------------------------------------------------------------------------

SA_outputs <- read.table("C:/Users/Maxime/Documents/Git/Flat_fish_2021/01_data/stock_assessment_outputs/outputs.csv", header=T, sep=",")
SA_outputs <- as_tibble(SA_outputs)
colnames(SA_outputs)[1] <- c("Years")
colnames(SA_outputs)[2] <- c("q")

SA_outputs <- SA_outputs %>% dplyr::filter(Years %in% c(2001:2019)) %>% mutate(Residuals =(Obs_surv-PredSurvey)/ Std_obs_survey)

Residuals <- SA_outputs %>% dplyr::select(Years, Residuals) %>% mutate(Type="Residuals")
Residuals <- Residuals %>% 
  #mutate(Index = (Residuals-mean(Residuals))/sd(Residuals)) %>% dplyr::select(Years,Index,Type)
  mutate(Index = Residuals) %>% dplyr::select(Years,Index,Type)
# ------------------------------------------------

O_index_int <- O_uvt %>% dplyr::filter(Seasons =="Int.")
O_index <- (O_index_int$Overlap - mean(O_index_int$Overlap))/sd(O_index_int$Overlap) 
O_index <- as.data.frame(O_index)

O_index_int <- cbind(O_index,O_index_int$Years)
colnames(O_index_int) <- c("Index", "Years")
O_index_int <- as_tibble(O_index_int) %>% mutate(Type = "Overlap")

Index_comp <- rbind(O_index_int,Residuals)

# no smooth
#----------------
pcompare <- ggplot(Index_comp, mapping=aes(x=Years, y=Index,color=Type),size=2) + geom_point(size=2)+ geom_line(size=1.5)+
  #geom_smooth(aes(color=Type,fill=Type))+
  #geom_point(Index_comp, mapping=aes(x=Years, y=IA_index, color=Env.Cond),size=4)+
  scale_color_manual(values=c( "purple","yellow3"))
p_comp_enc_res <- pcompare +  
  geom_rect(aes(xmin = 2000, xmax = 2001, ymin = -Inf, ymax = Inf),
            fill = "lightblue",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2001, xmax = 2006, ymin = -Inf, ymax = Inf),
            fill = "pink",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2006, xmax = 2014, ymin = -Inf, ymax = Inf),
            fill = "lightblue",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2014, xmax = 2017, ymin = -Inf, ymax = Inf),
            fill = "pink",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2017, xmax = 2018, ymin = -Inf, ymax = Inf),
            fill = "lightblue",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2018, xmax = 2019, ymin = -Inf, ymax = Inf),
            fill = "pink",col=NA, alpha = 0.01) + ggtitle('c)')

ggsave(paste0(work_dir,'overlap_vs_residuals.png'), p_comp_enc)



# smooth
#--------------
pcompare_s <- ggplot(Index_comp, mapping=aes(x=Years, y=Index,color=Type),size=2) + geom_point(size=2)+
  geom_smooth(aes(color=Type),se=FALSE,size=2)+
  scale_color_manual(values=c( "purple","yellow3"))+
  scale_fill_manual(values = c("purple", "yellow3"))  
pcompare_s


pcompare_s_enc <- pcompare_s +  
  geom_rect(aes(xmin = 2001, xmax = 2006, ymin = -Inf, ymax = Inf),
            fill = "pink",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2006, xmax = 2014, ymin = -Inf, ymax = Inf),
            fill = "lightblue",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2014, xmax = 2017, ymin = -Inf, ymax = Inf),
            fill = "pink",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2017, xmax = 2018, ymin = -Inf, ymax = Inf),
            fill = "lightblue",col=NA, alpha = 0.01)+
  geom_rect(aes(xmin = 2018, xmax = 2019, ymin = -Inf, ymax = Inf),
            fill = "pink",col=NA, alpha = 0.01) + ggtitle('c)')



pseason


p <- pseason  / (p_comp_enc + p_comp_enc_res)
ggsave(paste0(work_dir,'Overlap_catch_residuals.png'), p)
getwd()



# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

load(paste0("C:/Users/Maxime/Documents/Git/Flat_fish_2021/04_outputs/user_grid_outEBS/sex/female/no2000/O_uvt_female.RData"))

O <- bind_rows(O_uvt_male,O_uvt_female)


p <- ggplot() + 
  geom_line(O, mapping=aes(x=Years, y=Overlap_m,color=Env.Cond,size=var),alpha=0.2)+
  ylab("Overlap")+
  geom_line(O, mapping=aes(x=Years, y=Overlap),size=1)+
  geom_point(O, mapping=aes(x=Years, y=Overlap, color=Env.Cond),size=2)+ guides(size = "none")+
  scale_color_manual(values=c( "#1E90FF","#FF0000"))

pseason <- p + facet_grid( sex~ Seasons)
pseason

x <- O %>% filter(Seasons =="Int.")

(0.6634530 -0.6169117 )/0.6169117
(0.7823043-0.8389110)/0.7823043


unique(x$Overlap_m)
ggsave(paste0('C:/Users/Maxime/Documents/Git/Flat_fish_2021/04_outputs/user_grid_outEBS/sex/Overlap_malefemale.png'), pseason)
