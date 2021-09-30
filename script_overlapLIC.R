###########################################################################
# CALCULATE OVERLAP -----------------------------------------------------------
###########################################################################
rm(list=ls())

# Packages ----------------------------------------------------------------
library(reshape2)
library(dplyr)
library(ggplot2)


# -------------------------------------------------------------------------
# OVERLAP BETWEEN FISHERIES AND SURVEY ------------------------------------
# -------------------------------------------------------------------------

# Time = Year --------------------------------------------------------
# -------------------------------------------------------------------------

# - fisheries CPUE --------------------------------------------------------
load("04_outputs/EBS_grid/EnvCxSeason/fisheries/c(3,3,3,3,2,2)/fit.RData")
fit_f <- fit
dim(fit_f$Report$D_gct)

cold <- c(6,7,8,9,10,11,12,13,17)
cold_year <- c(2006,2007,2008,2009,2010,2011,2012,2013,2017)

warm <- c(1,2,3,4,5,14,15,16,18)
warm_year <- c(2001,2002,2003,2004,2005,2014,2015,2016,2018)

n_knot <- fit_f$spatial_list$n_x 
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
      D_cold_gut[g,u,t] <- fit_f$Report$D_gct[g,1,cold[t]] * exp(fit_f$Report$Phi1_sk[g,cold_season[u]]) 
    }
  }
}

for ( t in 1:length(warm)){
  for ( u in 1 :length(warm_season)){
    for ( g in 1:n_knot){
      D_warm_gut[g,u,t] <- fit_f$Report$D_gct[g,1,warm[t]] * exp(fit_f$Report$Phi1_sk[g,warm_season[u]])  
    }
  }
}


# Dplyr format ------------------------------------------------------------
D_warm_gut <- melt(D_warm_gut)
unique(D_warm_gut$Var3)
colnames(D_warm_gut) <- c("knot","season","year","value")
temp_warm<- cbind(seq(1,length(warm_year)),warm_year)
colnames(temp_warm) <- c("year","Year")
D_warm_gut <- left_join(D_warm_gut,as_tibble(temp_warm))
D_warm_gut <- D_warm_gut %>% mutate(EnvCond = "Warm")



D_cold_gut <- melt(D_cold_gut)
colnames(D_cold_gut) <- c("knot","season","year","value")
temp_cold<- cbind(seq(1,length(cold_year)),cold_year)
colnames(temp_cold) <- c("year","Year")
D_cold_gut <- left_join(D_cold_gut,as_tibble(temp_cold))
D_cold_gut <- D_cold_gut %>% mutate(EnvCond = "Cold")



D_f <- bind_rows(D_warm_gut,D_cold_gut)
unique(D_f$Year)


# - survey CPUE -----------------------------------------------------------
load("04_outputs/EBS_grid/EnvCxSeason/survey/fit.RData")
fit_s <- fit
dim(fit_s$Report$D_gct)

D_s <- melt(fit_s$Report$D_gct)
colnames(D_s) <- c("knot","cat","year", "value_s")
D_s <- D_s %>% mutate(Year= year+1981) %>% dplyr::select(-year)


# - Join Survey and fisheries ---------------------------------------------

D <- left_join(D_f ,D_s)


# LIC ---------------------------------------------------------------------

D_temp <- D %>% group_by(Year,EnvCond,season) %>% summarise(LIC=sum(value*value_s)/sqrt(sum(value_s*value_s)*sum(value*value)))
D_temp$Season <- as.factor(D_temp$season)

D_temp <- D_temp %>% mutate(Seasons=recode(Season, "1" ="Ealry","2" ="Intermediate", "3"="Late"))
p <- ggplot(D_temp)+
  geom_line(aes(x=Year,y=LIC,color=Seasons),size=2)+geom_point(aes(x=Year, y=LIC,size=EnvCond))+theme_bw()
p 


png(paste('LIC_seasons','.png',sep=''), height =14 , width = 20, units = 'cm', res=500)#, res=500)
p
dev.off()


