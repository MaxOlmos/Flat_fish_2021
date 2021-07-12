
###########################################################################
########           Script CPUE VAST 
###########################################################################

# -------------------------------------------------------------------------
# GOAL :  -----------------------------------------------------------------
# Using CPUE fisheries data to model change in phenology (=timing of migration
# from wintering areas = offshore to spawning area= inshore)

# How ? Use VAST to model seasonal changes in density implicitly
## - week and cold pool effect are used as catchability covariates
## - Cold pool effect is here represented by a discrete variable = cold or warm years

# We want here to model the spatially varying effect of (i) week effect, (ii) cold pool effect and (iii) interaction of 
# cold pool X week on probability of encounter

# -------------------------------------------------------------------------
rm(list=ls())
getwd()
# Load package ------------------------------------------------------------
# Install FishStatsUtils from CRAN

library(VAST)
library(TMB)
library(dplyr)
library(rgdal)
library(raster)
library(pander)
library(splines)
library(FishStatsUtils)
library(devtools)
#devtools::install_github("james-thorson/FishData")
library(FishData)
#Install package
#install_github("james-thorson/VAST", INSTALL_opts="--no-staged-install")
# Load package
library(VAST)

getwd()
DateFile <- "03_model/EBS_grid"
#setwd(DateFile)
outputs_file <- "04_outputs/EBS_grid/EnvCxSeason/survey" 


# PREPARE DATA ----------------------------------------------------------
# -------------------------------------------------------------------------
Database_EBS = FishData::download_catch_rates( survey="EBSBTS", species_set="Limanda aspera" )
head(Database_EBS)

Database_EBS %>% group_by(Year) %>% summarise(sum=sum(Wt))

# ^ Check NA --------------------------------------------------------------
Database_EBS %>% filter(is.na(Long))

# ^ Prepare Data for VAST -------------------------------------------------
Data_Geostat <- Database_EBS

colnames(Data_Geostat)<-c("Sci", "Year","TowID" , "Lat", "Lon","Wt", "AreaSwept_ha") 

# Spatial Settings --------------------------------------------------------
library(sp)

# ^ Number of knot/stations -----------------------------------------------
n_x = c(30, 50, 75, 100, 150, 200, 300)[4]

# ^ Output from Calc_Kmeans with knots for a triangulated mesh ------------
# ^ Calc_Kmeans determines the location for a set of knots for ------------
# ^ approximating spatial variation ---------------------------------------
# ^ n_x: the number of knots to select ------------------------------------
# ^ nstart the number of times that the k-means algorithm is run while  ---
# ^ searching for the best solution (default=100) -------------------------
# ^ iter.max the number of iterations used per k-means algorithm ----------
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )

strata.limits <- data.frame('STRATA'="All_areas")
#strata.limits <- data.frame('STRATA'=list(All_areas = 1:1e+05))
Region_EBS = c("Eastern_Bering_Sea")
Extrapolation_List_EBS  = FishStatsUtils::make_extrapolation_info(Region=Region_EBS, strata.limits=strata.limits)


# - Settings ----------------------------------------------------------------

settings = make_settings( n_x = 100,
                          Region = "Eastern_Bering_Sea",
                          fine_scale=FALSE,
                          purpose = "index2",
                          bias.correct = FALSE,
                          ObsModel = c(2,0),
                          strata.limits=strata.limits
                          )

settings$strata.limits

# - Switch off estimation of random effect for spatial effect  ---------------
# and spatio-temporal effect ----------------------------------------------
# -------------------------------------------------------------------------

settings$FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1) 

# * * structure on parameters among years --------------------------------------------------------------------
# strucutre for intercetps ie spatio temporal varition across time, uisng input
settings$RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)

# * * Number of oversdispersiom factors= Number of cathcability factors---------------------------------------
settings$OverdispersionConfig = c("Eta1"=0, "Eta2"=0)


# RUN MODEL ---------------------------------------------------------------
# -------------------------------------------------------------------------

# -- Covariates -----------------------------------------------------------
settings$Options["report_additional_variables"]=FALSE

check_fit(par)

# -- run model
# -------------------------------------------------------------------------
Data_Geostat <- Data_Geostat %>% mutate (AreaSwept_km2=AreaSwept_ha/100)
Data_Geostat <- Data_Geostat %>% mutate (Catch_KG = 1000*Wt) 
Data_Geostat <- Data_Geostat %>% filter(Year > 2000) %>% filter(Year < 2019 )

# RUN MODEL ---------------------------------------------------------------
names(Data_Geostat)
workdir <- paste0(getwd(),"/",outputs_file)


fit = fit_model( "settings"=settings, 
                 "Lat_i"=Data_Geostat[,'Lat'], 
                 "Lon_i"=Data_Geostat[,'Lon'],
                 "t_i"=Data_Geostat[,'Year'],
                 "c_i"=rep(0,nrow(Data_Geostat)),
                 "b_i"=Data_Geostat[,'Catch_KG'],
                 "a_i"=Data_Geostat[,'AreaSwept_km2'],
                 working_dir=workdir)





save(fit, file=paste0(workdir,"/fit.RData"))
load(paste0(workdir,"/fit.RData"))


# -- plot ------------------------------------------------------------------
# -------------------------------------------------------------------------
load(paste0("fit.RData"))
plot( fit,plot_set=c(3,6,7,15,16),
      Yrange=c(NA,NA),working_dir =paste0(workdir,"/"))


plot( fit,
      Yrange=c(NA,NA))


dharmaRes = summary(fit, what="residuals", working_dir=workdir)

file.info()

getwd()
# Test significativity ----------------------------------------------------
# -------------------------------------------------------------------------
# Re-run model with getJointPrecision=TRUE --------------------------------
fit_sign = fit_model( settings = settings,
                      Lat_i = CPUE_catchability[,'lat'],
                      Lon_i = CPUE_catchability[,'long'],
                      t_i = CPUE_catchability[,'year'],
                      b_i = CPUE_catchability[,'CPUE.mean'],
                      a_i = rep(1, nrow(CPUE_catchability)),
                      Map=Map,
                      Q1config_k =c(3,3,3,3,2,2), 
                      Q_ik=Q1_ik,
                      input_grid=user_region,
                      getJointPrecision=TRUE,
                      working_dir=workdir )

save(fit_sign, file=paste0(workdir,"/fit_sign.RData"))

load(paste0(workdir,"/fit_sign.RData"))


plot( fit_sign,
      Yrange=c(NA,NA),working_dir =paste0(workdir,"/"))


# Re sample  --------------------------------------------------------------
sample_true = sample_variable( Sdreport=fit_sign$parameter_estimates$SD,
                               Obj=fit_sign$tmb_list$Obj, variable_name="Phi1_sk" )
sample_false = sample_variable( Sdreport=fit_sign$parameter_estimates$SD,#n_samples = 30,
                                Obj=fit_sign$tmb_list$Obj, variable_name="Phi1_sk",sample_fixed=FALSE )

# TRUE
var <- apply(sample_true,c(1,2),var)
# FALSE-> work better with false
var <- apply(sample_false,c(1,2),var)
sd <- sqrt(var)
Wt <- fit$Report$Phi1_sk/sd # t value from t test in summary

# calcul pvalue : same as pvalue in summary() test
pvalue = pchisq(Wt^2, df = 1, lower = F) #
pvalue2 =pnorm(-abs(Wt), 0, 1, lower = T) + pnorm(abs(Wt), 0, 1, lower = F) #
round(pvalue-pvalue2, 5)

pvaluenNA <- ifelse(pvalue>0.05, pvalue==NA, pvalue)


# -------------------------------------------------------------------------

# I did not find the function to plot the spatially response of --------
# covariates so I made the plot by hand here
# -------------------------------------------------------------------------
getwd()
MapDetails_List <- FishStatsUtils:: make_map_info( "Region"=settings$Region,spatial_list = fit$spatial_list, "NN_Extrap"=fit$spatial_list$PolygonList$NN_Extrap,                                                     
                                                   "Extrapolation_List"=fit$extrapolation_list ,fine_scale = FALSE)

shape_path <- "C:/Users/Maxime/Documents/Git/Flat_fish/Pheno-flatfish/ADFGv2/shapefile/"
coast_shapefile <- paste(shape_path, "ne_50m_land.shp", sep="")
ocean <- readOGR(coast_shapefile)

head(Q1_ik)



# Plot covariates effects -------------------------------------------------

phi <- fit$Report$Phi1_sk


bathy.dat<-read.table('C:/Users/Maxime/Documents/Git/Flat_fish_2021/01_data/bathyEBS/BeringDepth.txt',sep='') 
names(bathy.dat)<-c('lon','lat','depth')

bathy.dat$depth[bathy.dat$depth>0]<-NA#Avoid points above water
head(bathy.dat)
bathy.mat<- matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]


png(paste(workdir,"/Covariates_effect",'.png',sep=''), height = 8, width = 12, units = 'in', res=600)
par(mfrow=c(2,3))
#par(mar=c(4.5, 4.5, 0.5, 0.5))
name_effect <- c("ColdEarly","WarmEarly","ColdInt","WarmInt","ColdLate","WarmLate")#"ColdLate","WarmLate")
Xplot <- pvaluenNA#fit$Report$Phi1_sk[,]
CP<- c("cold","warm","cold","warm","cold","warm")
Season <- c("early","early","int","int","late","late")
Season <- as.character(Season)

for (k in c(1,3,5,2,4,6)){
  
  x <- Season[k]
  y<- CP[k]
  data_obs <- as_tibble(CPUE_catchability) %>% dplyr:: filter(Season == x, CP==y)
  
  x <- as.data.frame(MapDetails_List$PlotDF)
  PlotDF1 = subset(x,Include==TRUE)
  
  unique(x$Include)
  PlotDF1 <- as.data.frame(PlotDF1)
  coords2 = cbind(PlotDF1$Lon,PlotDF1$Lat)
  coords2 = coords2[coords2[,1]<0,]
  P2 = SpatialPoints(coords2)
  P3 = SpatialPoints(coords2)
  
  rast <- raster(ncol=50,nrow=50)
  extent(rast) <- extent(coords2)
  col=colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                         "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))#colorRampPalette( c('blue','red', 'yellow'), bias=1)
  
  breakpoints  <-  seq(min( fit$Report$Phi1_sk[PlotDF1$x2i,]),round(max(fit$Report$Phi1_sk[PlotDF1$x2i,]),digits=1),(round(max(fit$Report$Phi1_sk[PlotDF1$x2i,]),digits=1)-min(fit$Report$Phi1_sk[PlotDF1$x2i,]))/70)
  range(fit$Report$Phi1_sk)
  
  P2$data =  (phi[,k])[PlotDF1$x2i]
  rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
  
  
  
  plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(54,62),main=name_effect[k])
  
  image(rast.temp,col=col(length(breakpoints)-1),axes=TRUE,breaks=breakpoints,
        add=T,xlim=c(-180,-158),ylim=c(54,62))
  contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(200,100,50),labcex=0.4,col='black',add=T)
  points(data_obs$long,data_obs$lat, pch=3)
  
  box()
  legend_x=c(0.1,0.2)
  legend_y=c(0.05,0.3)
  cex.legend=0.4
  xl = (1-legend_x[1])*par('usr')[1] + (legend_x[1])*par('usr')[2]
  xr = (1-legend_x[2])*par('usr')[1] + (legend_x[2])*par('usr')[2]
  yb = (1-legend_y[1])*par('usr')[3] + (legend_y[1])*par('usr')[4]
  yt = (1-legend_y[2])*par('usr')[3] + (legend_y[2])*par('usr')[4]
  if( diff(legend_y) > diff(legend_x) ){
    align = c("lt","rb")[2]
    gradient = c("x","y")[2]
  }else{
    align = c("lt","rb")[1]
    gradient = c("x","y")[1]
  }  
  
  plotrix::color.legend(xl=xl, yb=yb, xr=xr, yt=yt, legend=round(seq(min(breakpoints),max(breakpoints),length=4),2),
                        rect.col=col(length(breakpoints)-1), cex=0.6, align=align, gradient=gradient)
  
  
  
}

dev.off()

# Highligting significant phi effects -------------------------------------
# -------------------------------------------------------------------------

for (i in 1:dim(phi)[1]){
  for (j in 1:dim(phi)[2]){
    if(pvalue[i,j]>0.05){phi[i,j]=NA}else{phi[i,j]==phi[i,j]}
    
  }
}


png(paste(workdir,"/Covariates_significanteffect",'.png',sep=''), height = 8, width = 12, units = 'in', res=600)
par(mfrow=c(2,3))
#par(mar=c(4.5, 4.5, 0.5, 0.5))
name_effect <- c("ColdEarly","WarmEarly","ColdInt","WarmInt","ColdLate","WarmLate")#"ColdLate","WarmLate")
Xplot <- pvaluenNA#fit$Report$Phi1_sk[,]
CP<- c("cold","warm","cold","warm","cold","warm")
Season <- c("early","early","int","int","late","late")
Season <- as.character(Season)

for (k in c(1,3,5,2,4,6)){
  
  x <- Season[k]
  y<- CP[k]
  data_obs <- as_tibble(CPUE_catchability) %>% dplyr:: filter(Season == x, CP==y)
  
  x <- as.data.frame(MapDetails_List$PlotDF)
  PlotDF1 = subset(x,Include==TRUE)
  
  unique(x$Include)
  PlotDF1 <- as.data.frame(PlotDF1)
  coords2 = cbind(PlotDF1$Lon,PlotDF1$Lat)
  coords2 = coords2[coords2[,1]<0,]
  P2 = SpatialPoints(coords2)
  P3 = SpatialPoints(coords2)
  
  rast <- raster(ncol=50,nrow=50)
  extent(rast) <- extent(coords2)
  col=colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                         "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))#colorRampPalette( c('blue','red', 'yellow'), bias=1)
  
  breakpoints  <-  seq(min( fit$Report$Phi1_sk[PlotDF1$x2i,]),round(max(fit$Report$Phi1_sk[PlotDF1$x2i,]),digits=1),(round(max(fit$Report$Phi1_sk[PlotDF1$x2i,]),digits=1)-min(fit$Report$Phi1_sk[PlotDF1$x2i,]))/70)
  range(fit$Report$Phi1_sk)
  
  P2$data =  (phi[,k])[PlotDF1$x2i]
  rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
  
  P3$data =  (pvalue[,k])[PlotDF1$x2i]
  rast.temp3 <- rasterize(P3, rast, P3$data, fun = mean)
  
  
  
  plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(54,62),main=name_effect[k])
  image(rast.temp3,col="lightgrey",axes=TRUE,
        add=T,xlim=c(-180,-158),ylim=c(54,62))
  image(rast.temp,col=col(length(breakpoints)-1),axes=TRUE,breaks=breakpoints,
        add=T,xlim=c(-180,-158),ylim=c(54,62))
  contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(200,100,50),labcex=0.4,col='black',add=T)
  points(data_obs$long,data_obs$lat, pch=3)
  
  box()
  legend_x=c(0.1,0.2)
  legend_y=c(0.05,0.3)
  cex.legend=0.4
  xl = (1-legend_x[1])*par('usr')[1] + (legend_x[1])*par('usr')[2]
  xr = (1-legend_x[2])*par('usr')[1] + (legend_x[2])*par('usr')[2]
  yb = (1-legend_y[1])*par('usr')[3] + (legend_y[1])*par('usr')[4]
  yt = (1-legend_y[2])*par('usr')[3] + (legend_y[2])*par('usr')[4]
  if( diff(legend_y) > diff(legend_x) ){
    align = c("lt","rb")[2]
    gradient = c("x","y")[2]
  }else{
    align = c("lt","rb")[1]
    gradient = c("x","y")[1]
  }  
  
  plotrix::color.legend(xl=xl, yb=yb, xr=xr, yt=yt, legend=round(seq(min(breakpoints),max(breakpoints),length=4),2),
                        rect.col=col(length(breakpoints)-1), cex=0.6, align=align, gradient=gradient)
  
  
  
}

dev.off()
