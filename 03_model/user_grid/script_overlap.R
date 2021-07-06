###########################################################################
# CALCULATE OVERLAP -----------------------------------------------------------
###########################################################################
rm(list=ls())

# Packages ----------------------------------------------------------------
library(reshape2)
library(dplyr)
library(ggplot2)

# Load Data ---------------------------------------------------------------
DateFile <- "03_model/user_grid"
outputs_file <- "04_outputs/user_grid/EnvCxSeason/c(3,3,3,3,2,2)_grid1" 

workdir <- paste0(getwd(),"/",outputs_file)

load(paste0(workdir,"/fit.RData"))


# Calculating COG for cold/warm Years per Seasons -------------------------
# - Calculating Densitity for cold and warm year and season ---------------
dim((fit$Report$D_gct))
dim(fit$Report$Phi1_sk)

cold <- c(6,7,8,9,10,11,12,13,17)
cold_year <- c(2006,2007,2008,2009,2010,2011,2012,2013,2017)

warm <- c(1,2,3,4,5,14,15,16,18)
warm_year <- c(2001,2002,2003,2004,2005,2014,2015,2016,2018)

n_knot <- fit$spatial_list$n_x 
name_effect <- c("ColdEarly","WarmEarly","ColdInt","WarmInt","ColdLate","WarmLate")
season <- c("Early","Int.","Late")
cold_season <- c(1,3,5)
warm_season <- c(2,4,6)
D_cold_gut <- array(NA,dim=c(n_knot,length(cold_season),length(cold)))
dim(D_cold_gut)
D_warm_gut <- array(NA,dim=c(n_knot,length(warm_season),length(warm)))

for ( t in 1:length(cold)){
  for ( u in 1 :length(cold_season)){
    for ( g in 1:n_knot){
      D_cold_gut[g,u,t] <- fit$Report$D_gct[g,1,cold[t]] * exp(fit$Report$Phi1_sk[g,cold_season[u]]) 
    }
  }
}

for ( t in 1:length(warm)){
  for ( u in 1 :length(warm_season)){
    for ( g in 1:n_knot){
      D_warm_gut[g,u,t] <- fit$Report$D_gct[g,1,warm[t]] * exp(fit$Report$Phi1_sk[g,warm_season[u]])  
    }
  }
}

D_cold_gu_plot <- D_cold_gut
D_warm_gut_plot <- D_warm_gut

# - Average accross all years ---------------------------------------------
D_cold_gut_m <- apply(D_cold_gut,c(1,2),mean)
D_warm_gut_m <- apply(D_warm_gut,c(1,2),mean)

# Dplyr format ------------------------------------------------------------
D_warm_gut <- melt(D_warm_gut)
colnames(D_warm_gut) <- c("knot","season","year","value")

D_cold_gut <- melt(D_cold_gut)
colnames(D_cold_gut) <- c("knot","season","year","value")

loc <- cbind(fit$spatial_list$loc_x, seq(1,n_knot,1))
colnames(loc) <-c("E_km", "N_km", "knot")
loc <- as.data.frame(loc)

fit$spatial_list$loc_i

D_warm_gut <- left_join(D_warm_gut, loc, by= "knot")
D_cold_gut <- left_join(D_cold_gut, loc, by= "knot")

# - Calculate abundances --------------------------------------------------
A_g <- fit$spatial_list$a_gl*fit$spatial_list$a_gl*1000
A_g <- cbind(A_g, seq(1,n_knot,1))
colnames(A_g) <- c("area_km2","knot")
A_g <- as.data.frame(A_g)

D_warm_gut <- left_join(D_warm_gut,A_g) 
D_cold_gut <- left_join(D_cold_gut,A_g) 

D_warm_gut <-  D_warm_gut %>% mutate(Ab = value*area_km2)
D_cold_gut <-  D_cold_gut %>% mutate(Ab = value*area_km2)