
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


# PREPARE DATA ----------------------------------------------------------
# -------------------------------------------------------------------------
load(file="20022_transformed_data/observer/CPUE_catchability_adfg.RData")
CPUE_catchability
Sample_size <- CPUE_catchability %>% dplyr::group_by(year, Season, CP) %>% dplyr::summarise(sample=n())
ggplot()+ geom_line(data=Sample_size, mapping = aes(x=year,y=sample, color = Season))


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
settings$FieldConfig[1,1]<-0
settings$FieldConfig[2,1]<-0

#if 2 kappas
#settings$FieldConfig[1,2]<-0
#settings$FieldConfig[2,2]<-0
#settings$FieldConfig[2,1] <- 0



# RUN MODEL ---------------------------------------------------------------
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# - Model : SeasonXCP covariate -----------------------------------------------------
# -------------------------------------------------------------------------
getwd()
outputs_file <- "2022/04_outputs/user_grid_outEBS/no_sex/M3s" 

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
??fit_model
fit = fit_model( settings = settings,
                 Lat_i = CPUE_catchability[,'lat'],
                 Lon_i = CPUE_catchability[,'long'],
                 t_i = CPUE_catchability[,'year'],
                 b_i = CPUE_catchability[,'CPUE.mean'],
                 a_i = rep(1, nrow(CPUE_catchability)),
                 Q2config_k = c(3,3,3,3,2,2), 
                 Q_ik=Q2_ik,build_model = FALSE,
                 input_grid=user_region,working_dir=workdir)
fit$data_list$n_g

# plot knots
knot_loc2 <- as.data.frame(fit$spatial_list$latlon_g)
plot(knot_loc2$Lon,knot_loc2$Lat, col="red")


# Modify Map
Map = fit$tmb_list$Map
Map$lambda1_k = rep(factor(NA),6)
#Map$lambda1_k = factor(c(1,1,2,2,3,3))
Map$log_sigmaPhi2_k <- factor(rep(as.numeric(Map$log_sigmaPhi2_k[1]),6))
getwd()

fit = fit_model( settings = settings,
                 Lat_i = CPUE_catchability[,'lat'],
                 Lon_i = CPUE_catchability[,'long'],
                 t_i = CPUE_catchability[,'year'],
                 b_i = CPUE_catchability[,'CPUE.mean'],
                 a_i = rep(1, nrow(CPUE_catchability)),
                 Map=Map,
                 Q2config_k =c(3,3,3,3,2,2), 
                 Q_ik=Q2_ik,
                 input_grid=user_region,
                 #getJointPrecision=TRUE,
                 working_dir=workdir )

save(fit, file=paste0(workdir,"/fit.RData"))
#load(paste0(workdir,"/fit.RData"))
getwd()
fit$parameter_estimates$par
# -- plot ------------------------------------------------------------------
# -------------------------------------------------------------------------

plot( fit,plot_set=c(3,6,7,15,16,20),
      Yrange=c(NA,NA),working_dir =paste0(workdir,"/"))



plot( fit,
      Yrange=c(NA,NA))


dharmaRes = summary(fit, what="residuals", working_dir=workdir)


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
                      Q2config_k =c(3,3,3,3,2,2), 
                      Q_ik=Q2_ik,
                      input_grid=user_region,
                      getJointPrecision=TRUE,
                      working_dir=workdir )

save(fit_sign, file=paste0(workdir,"/fit_sign.RData"))

#load(paste0(workdir,"/fit_sign.RData"))


plot( fit_sign,
      Yrange=c(NA,NA),working_dir= workdir)


# Re sample  --------------------------------------------------------------
sample_true = sample_variable( Sdreport=fit_sign$parameter_estimates$SD,
                               Obj=fit_sign$tmb_list$Obj, variable_name="Phi2_gk",n_samples = 250 )
sample_false = sample_variable( Sdreport=fit_sign$parameter_estimates$SD,n_samples = 250,
                                Obj=fit_sign$tmb_list$Obj, variable_name="Phi2_gk",sample_fixed=FALSE )

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

