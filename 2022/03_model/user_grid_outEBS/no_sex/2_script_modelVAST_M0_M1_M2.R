
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
library(ggplot2)
library(sf)
library(ne)
# Install package
#install_github("james-thorson/VAST", INSTALL_opts="--no-staged-install")
# Load package
library(VAST)
library(viridis)
version= get_latest_version(package="VAST")

# PREPARE DATA ----------------------------------------------------------
# -------------------------------------------------------------------------
load(file="02_transformed_data/observer/CPUE_catchability_adfg.RData")
CPUE_catchability
Sample_size <- CPUE_catchability %>% dplyr::group_by(year, Season, CP) %>% dplyr::summarise(sample=n())
ggplot()+ geom_line(data=Sample_size, mapping = aes(x=year,y=sample, color = Season))

###############
map <- ggplot()+ 
  geom_point(data=CPUE_catchability,aes(x=long,y=lat,col=log(CPUE.mean)),size=1)+
  #geom_sf(data=Ab_map %>% filter(size_class==1),aes(col=log(A_pkt)))+
  scale_color_viridis() + labs(colour= "Observed CPUE (log)")+
  #   geom_sf(data=EBS,fill=NA, color="brown")+
#  geom_sf(data = sebs_layers$bathymetry,size=0.2) +
  #geom_sf(data=world, col=NA, fill="black")+
  #coord_sf(xlim = xlims, ylim = ylims)+
 # facet_grid(year ~ Season)+
  theme(axis.text = element_blank())
#theme_bw(base_size = 8)

ggsave((paste0(work_dir, "data2.png")),plot=map,
       width = 27,
       height = 30,
       units = "cm")
##################

user_region <- readRDS('02_transformed_data/user_grid_outEBS/user_region.rds')
plot(user_region$Lon,user_region$Lat)
#user_region$Area_km2 <- 1

# - Settings ----------------------------------------------------------------
settings = make_settings( n_x = 100,
                          Region = "User",
                          #fine_scale=FALSE,
                          purpose = "index2",
                          bias.correct = FALSE,
                          ObsModel = c(1,4)#,
                          #knot_method='grid'
                          #knot_method='samples'
)
settings$fine_scale
??make_settings
#fit$settings$knot_method

# - Switch off estimation of random effect for spatial effect  ---------------
# and spatio-temporal effect ----------------------------------------------
# -------------------------------------------------------------------------

## WIHTOUT COVARIATES
# because variance paramters are apporaching zero
# because variance paramters are apporaching zero
#settings$FieldConfig[1,2]<-0
#settings$FieldConfig[2,2]<-0
#settings$FieldConfig[2,1] <- 0

settings$FieldConfig[1,1] <- 0
settings$FieldConfig[2,1] <- 0

# RUN MODEL ---------------------------------------------------------------
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# - Model : SeasonXCP covariate -----------------------------------------------------
# -------------------------------------------------------------------------
getwd()
outputs_file <-  "2022/04_outputs/user_grid_outEBS/no_sex/M0" 

png(paste0(outputs_file,'/Interractions.png'), height = 6, width = 10, units = 'in', res=600)
par(mfrow=c(1,2))
with(CPUE_catchability,interaction.plot(Season, CP, log(CPUE.mean)))
with(CPUE_catchability,interaction.plot(CP, Season, log(CPUE.mean)))
dev.off()

# -- Covariates -----------------------------------------------------------
settings$Options["report_additional_variables"]=TRUE

# Levels
table(CPUE_catchability$CP)
CPUE_catchability <- CPUE_catchability %>% filter(year>2000)
# Build model matrix
Q1_formula = ~   CP:Season
Model_matrix1 = model.matrix( update.formula(Q1_formula, ~.), data=CPUE_catchability)
Columns_to_keep = which( attr(Model_matrix1,"assign") >0 )
coefficient_names_Q1 = attr(Model_matrix1,"dimnames")[[2]][Columns_to_keep]
Q2_ik = Model_matrix1[,Columns_to_keep,drop=FALSE]
head(Q1_ik)
apply(Q1_ik,2,sum)

# -- run model
# -------------------------------------------------------------------------

CPUE_catchability <- as.data.frame(CPUE_catchability)
workdir <- paste0(getwd(),"/",outputs_file,"/")

fit = fit_model( settings = settings,
                 Lat_i = CPUE_catchability[,'lat'],
                 Lon_i = CPUE_catchability[,'long'],
                 t_i = CPUE_catchability[,'year'],
                 b_i = CPUE_catchability[,'CPUE.mean'],
                 a_i = rep(1, nrow(CPUE_catchability)),
                 #Q1config_k = c(3,3,3,3,2,2), 
                 #Q_ik=Q1_ik,build_model = FALSE,
                 input_grid=user_region,working_dir=workdir,test_fit=FALSE)


save(fit, file=paste0(workdir,"/fit.RData"))
#load(paste0(workdir,"/fit.RData"))
getwd()


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
#M1
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------



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
version= get_latest_version(package="VAST")

# PREPARE DATA ----------------------------------------------------------
# -------------------------------------------------------------------------
load(file="2022/02_transformed_data/observer/CPUE_catchability_adfg.RData")

user_region <- readRDS('02_transformed_data/user_grid_outEBS/user_region.rds')
plot(user_region$Lon,user_region$Lat)
#user_region$Area_km2 <- 1

# - Settings ----------------------------------------------------------------
getwd()
outputs_file <-"2022/04_outputs/user_grid_outEBS/no_sex/M1" 
workdir <- paste0(getwd(),"/",outputs_file,"/")

settings = make_settings( n_x = 100,
                          Region = "User",
                          #fine_scale=FALSE,
                          purpose = "index2",
                          bias.correct = FALSE,
                          ObsModel = c(1,4)#,
                          #knot_method='grid'
                          #knot_method='samples'
)

#fit$settings$knot_method

# - Switch off estimation of random effect for spatial effect  ---------------
# and spatio-temporal effect ----------------------------------------------
# -------------------------------------------------------------------------

## WIHTOUT COVARIATES
# because variance paramters are apporaching zero
settings$FieldConfig[1,1] <- 0
settings$FieldConfig[2,1] <- 0


# -- Covariates -----------------------------------------------------------
settings$Options["report_additional_variables"]=TRUE

# Levels
table(CPUE_catchability$CP)
CPUE_catchability <- CPUE_catchability %>% filter(year>2000)
# Build model matrix
Q1_formula = ~   CP
Model_matrix1 = model.matrix( update.formula(Q1_formula, ~.), data=CPUE_catchability)
Columns_to_keep = which( attr(Model_matrix1,"assign") >0 )
coefficient_names_Q1 = attr(Model_matrix1,"dimnames")[[2]][Columns_to_keep]
Q2_ik = Model_matrix1[,Columns_to_keep,drop=FALSE]
head(Q2_ik)
CPCold = 1- Q2_ik[,1]
Q2_ik = cbind(Q2_ik,CPCold)

# -- run model
# -------------------------------------------------------------------------

CPUE_catchability <- as.data.frame(CPUE_catchability)

fit = fit_model( settings = settings,
                 Lat_i = CPUE_catchability[,'lat'],
                 Lon_i = CPUE_catchability[,'long'],
                 t_i = CPUE_catchability[,'year'],
                 b_i = CPUE_catchability[,'CPUE.mean'],
                 a_i = rep(1, nrow(CPUE_catchability)),
                 Q2config_k = c(2,2), 
                 Q_ik=Q2_ik,build_model = FALSE,
                 input_grid=user_region,working_dir=workdir)

# plot knots
knot_loc2 <- as.data.frame(fit$spatial_list$latlon_g)
plot(knot_loc2$Lon,knot_loc2$Lat, col="red")


# Modify Map
Map = fit$tmb_list$Map
Map$lambda1_k = rep(factor(NA),2)
#Map$lambda1_k = factor(c(1,1,2,2,3,3))
Map$log_sigmaPhi2_k <- factor(rep(as.numeric(Map$log_sigmaPhi2_k[1]),2))


fit = fit_model( settings = settings,
                 Lat_i = CPUE_catchability[,'lat'],
                 Lon_i = CPUE_catchability[,'long'],
                 t_i = CPUE_catchability[,'year'],
                 b_i = CPUE_catchability[,'CPUE.mean'],
                 a_i = rep(1, nrow(CPUE_catchability)),
                 Map=Map,
                 Q2config_k =c(2,2), 
                 Q_ik=Q2_ik,
                 input_grid=user_region,
                 #getJointPrecision=TRUE,
                 working_dir=workdir,test_fit=FALSE )

save(fit, file=paste0(workdir,"/fit.RData"))
#load(paste0(workdir,"/fit.RData"))
getwd()

settings$FieldConfig["Omega","Component_1"]

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
#M2
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------


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
version= get_latest_version(package="VAST")

getwd()
DateFile <- "03_model/user_grid"
#setwd(DateFile)

# PREPARE DATA ----------------------------------------------------------
# -------------------------------------------------------------------------
load(file="02_transformed_data/observer/CPUE_catchability_adfg.RData")
CPUE_catchability



user_region <- readRDS('02_transformed_data/user_grid_outEBS/user_region.rds')
plot(user_region$Lon,user_region$Lat)
#user_region$Area_km2 <- 1

# - Settings ----------------------------------------------------------------
settings = make_settings( n_x = 100,
                          Region = "User",
                          #fine_scale=FALSE,
                          purpose = "index2",
                          bias.correct = FALSE,
                          ObsModel = c(1,4)#,
                          #knot_method='grid'
                          #knot_method='samples'
)
settings$fine_scale
??make_settings
#fit$settings$knot_method

# - Switch off estimation of random effect for spatial effect  ---------------
# and spatio-temporal effect ----------------------------------------------
# -------------------------------------------------------------------------

## WIHTOUT COVARIATES
# because variance paramters are apporaching zero
settings$FieldConfig[1,1] <- 0
settings$FieldConfig[2,1] <- 0

# RUN MODEL ---------------------------------------------------------------
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# - Model : SeasonXCP covariate -----------------------------------------------------
# -------------------------------------------------------------------------
getwd()
outputs_file <- "2022/04_outputs/user_grid_outEBS/no_sex/M2" 

png(paste0(outputs_file,'/Interractions.png'), height = 6, width = 10, units = 'in', res=600)
par(mfrow=c(1,2))
with(CPUE_catchability,interaction.plot(Season, CP, log(CPUE.mean)))
with(CPUE_catchability,interaction.plot(CP, Season, log(CPUE.mean)))
dev.off()

# -- Covariates -----------------------------------------------------------
settings$Options["report_additional_variables"]=TRUE

# Levels
table(CPUE_catchability$CP)
CPUE_catchability <- CPUE_catchability %>% filter(year>2000)
# Build model matrix
Q1_formula = ~   Season
Model_matrix1 = model.matrix( update.formula(Q1_formula, ~.), data=CPUE_catchability)
Columns_to_keep = which( attr(Model_matrix1,"assign") >0 )
coefficient_names_Q1 = attr(Model_matrix1,"dimnames")[[2]][Columns_to_keep]
Q1_ik = Model_matrix1[,Columns_to_keep,drop=FALSE]

Seasonearly = 1- Q1_ik[,1]-Q1_ik[,2]
Q1_ik = cbind(Q1_ik,Seasonearly)
Q1_ik <- cbind(Q1_ik[,3],Q1_ik[,1],Q1_ik[,2])

# -- run model
# -------------------------------------------------------------------------

CPUE_catchability <- as.data.frame(CPUE_catchability)
workdir <- paste0(getwd(),"/",outputs_file,"/")

fit = fit_model( settings = settings,
                 Lat_i = CPUE_catchability[,'lat'],
                 Lon_i = CPUE_catchability[,'long'],
                 t_i = CPUE_catchability[,'year'],
                 b_i = CPUE_catchability[,'CPUE.mean'],
                 a_i = rep(1, nrow(CPUE_catchability)),
                 Q2config_k = c(3,3,2), 
                 Q_ik=Q1_ik,build_model = FALSE,
                 input_grid=user_region,working_dir=workdir)
fit$data_list$n_g

# plot knots
knot_loc2 <- as.data.frame(fit$spatial_list$latlon_g)
plot(knot_loc2$Lon,knot_loc2$Lat, col="red")


# Modify Map
Map = fit$tmb_list$Map
Map$lambda1_k = rep(factor(NA),3)
#Map$lambda1_k = factor(c(1,1,2,2,3,3))
Map$log_sigmaPhi2_k <- factor(rep(as.numeric(Map$log_sigmaPhi2_k[1]),3))


fit = fit_model( settings = settings,
                 Lat_i = CPUE_catchability[,'lat'],
                 Lon_i = CPUE_catchability[,'long'],
                 t_i = CPUE_catchability[,'year'],
                 b_i = CPUE_catchability[,'CPUE.mean'],
                 a_i = rep(1, nrow(CPUE_catchability)),
                 Map=Map,
                 Q2config_k =c(3,3,2), 
                 Q_ik=Q1_ik,
                 input_grid=user_region,
                 #getJointPrecision=TRUE,
                 working_dir=workdir,test_fit=FALSE )

save(fit, file=paste0(workdir,"/fit.RData"))
#load(paste0(workdir,"/fit.RData"))
getwd()


