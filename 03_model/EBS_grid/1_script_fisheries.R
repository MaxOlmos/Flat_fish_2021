
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
# Install package
#install_github("james-thorson/VAST", INSTALL_opts="--no-staged-install")
# Load package
library(VAST)

getwd()
DateFile <- "03_model/EBS_grid"
#setwd(DateFile)


# PREPARE DATA ----------------------------------------------------------
# -------------------------------------------------------------------------
load(file="02_transformed_data/observer/CPUE_catchability_adfg.RData")

# - Settings ----------------------------------------------------------------

settings = make_settings( n_x = 100,
                          Region = "Eastern_Bering_Sea",
                          fine_scale=FALSE,
                          purpose = "index2",
                          bias.correct = FALSE,
                          ObsModel = c(1,4) )

# - Switch off estimation of random effect for spatial effect  ---------------
# and spatio-temporal effect ----------------------------------------------
# -------------------------------------------------------------------------

## WIHTOUT COVARIATES
# because variance paramters are apporaching zero
settings$FieldConfig[1,2]<-0
settings$FieldConfig[2,1]<-0
#settings$FieldConfig[1,1]<-0


# RUN MODEL ---------------------------------------------------------------
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# - Model : SeasonXCP covariate -----------------------------------------------------
# -------------------------------------------------------------------------
getwd()
outputs_file <- "04_outputs/EBS_grid/EnvCxSeason/fisheries/c(3,3,3,3,2,2)" 

png(paste0(outputs_file,'/Interractions.png'), height = 6, width = 10, units = 'in', res=600)
par(mfrow=c(1,2))
with(CPUE_catchability,interaction.plot(Season, CP, log(CPUE.mean)))
with(CPUE_catchability,interaction.plot(CP, Season, log(CPUE.mean)))
dev.off()

# -- Covariates -----------------------------------------------------------
settings$Options["report_additional_variables"]=TRUE

# Levels
table(CPUE_catchability$CP)

# Build model matrix
Q1_formula = ~   CP:Season
Model_matrix1 = model.matrix( update.formula(Q1_formula, ~.), data=CPUE_catchability)
Columns_to_keep = which( attr(Model_matrix1,"assign") >0 )
coefficient_names_Q1 = attr(Model_matrix1,"dimnames")[[2]][Columns_to_keep]
Q1_ik = Model_matrix1[,Columns_to_keep,drop=FALSE]
head(Q1_ik)
apply(Q1_ik,2,sum)

# -- run model
# -------------------------------------------------------------------------

CPUE_catchability <- as.data.frame(CPUE_catchability)
workdir <- paste0(getwd(),"/",outputs_file)

fit = fit_model( settings = settings,
                 Lat_i = CPUE_catchability[,'lat'],
                 Lon_i = CPUE_catchability[,'long'],
                 t_i = CPUE_catchability[,'year'],
                 b_i = CPUE_catchability[,'CPUE.mean'],
                 a_i = rep(1, nrow(CPUE_catchability)),
                 Q1config_k = c(3,3,3,3,2,2), 
                 Q_ik=Q1_ik,build_model = FALSE,
                 working_dir=workdir)



# Modify Map
Map = fit$tmb_list$Map
Map$lambda2_k = rep(factor(NA),6)
#Map$lambda1_k = factor(c(1,1,2,2,3,3))
Map$log_sigmaPhi1_k <- factor(rep(as.numeric(Map$log_sigmaPhi1_k[1]),6))


fit = fit_model( settings = settings,
                 Lat_i = CPUE_catchability[,'lat'],
                 Lon_i = CPUE_catchability[,'long'],
                 t_i = CPUE_catchability[,'year'],
                 b_i = CPUE_catchability[,'CPUE.mean'],
                 a_i = rep(1, nrow(CPUE_catchability)),
                 Map=Map,
                 Q1config_k =c(3,3,3,3,2,2), 
                 Q_ik=Q1_ik,
                 #getJointPrecision=TRUE,
                 working_dir=workdir )

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
