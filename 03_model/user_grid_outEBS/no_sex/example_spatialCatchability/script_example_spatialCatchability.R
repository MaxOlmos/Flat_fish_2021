
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
# Load package
library(VAST)
version= get_latest_version(package="VAST")
getwd()

# LOAD DATA ----------------------------------------------------------
# -------------------------------------------------------------------------
load(file="CPUE_catchability_adfg.RData")

user_region <- readRDS('user_region.rds')
plot(user_region$Lon,user_region$Lat)

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

# - Switch off estimation of random effect for spatial effect  ---------------
# and spatio-temporal effect ----------------------------------------------
# -------------------------------------------------------------------------

## WIHTOUT COVARIATES
# because we are dealing with only positive CPUE
settings$FieldConfig[1,2]<-0
settings$FieldConfig[2,1]<-0


# RUN MODEL ---------------------------------------------------------------
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# - Model : SeasonXCP covariate -----------------------------------------------------
# -------------------------------------------------------------------------
getwd()
outputs_file <- ""## choose your file 

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
                 input_grid=user_region,working_dir=workdir)


# Modify Map
Map = fit$tmb_list$Map
Map$lambda2_k = rep(factor(NA),6)
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
                 input_grid=user_region,
                 #getJointPrecision=TRUE,
                 working_dir=workdir )

save(fit, file=paste0(workdir,"/fit.RData"))

# -- plot ------------------------------------------------------------------
# -------------------------------------------------------------------------

plot( fit,plot_set=c(3,6,7,15,16,20),
      Yrange=c(NA,NA),working_dir =paste0(workdir,"/"))



# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
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


plot( fit_sign,
      Yrange=c(NA,NA),working_dir= workdir)


# Re sample  --------------------------------------------------------------
sample_true = sample_variable( Sdreport=fit_sign$parameter_estimates$SD,
                               Obj=fit_sign$tmb_list$Obj, variable_name="Phi1_sk" )
sample_false = sample_variable( Sdreport=fit_sign$parameter_estimates$SD,n_samples = 250,
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

# test if it is the same to calcul pvalue from chisq or pnorm
round(pvalue-pvalue2, 5)

# prepare pvalue for plotting
pvaluenNA <- ifelse(pvalue>0.05, pvalue==NA, pvalue)


