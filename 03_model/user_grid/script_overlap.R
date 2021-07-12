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


# -------------------------------------------------------------------------
# OVERLAP WIHTIN FISHIERS CPUE BETWEEN SEASONS ----------------------------
# -------------------------------------------------------------------------


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
D_warm_gut_m <- melt(D_warm_gut_m)
colnames(D_warm_gut_m) <- c("knot","season","Warm")

D_cold_gut_m <- melt(D_cold_gut_m)
colnames(D_cold_gut_m) <- c("knot","season","Cold")

D_gut <- left_join(D_warm_gut_m,D_cold_gut_m)

# Calculate overlap -------------------------------------------------------

O_temp <- D_gut %>% mutate(numerator_0 =Cold*Warm ) 

O <- O_temp %>% group_by(season) %>% summarise(O = (2*sum(numerator_0))/(sum(Cold*Cold)+sum(Warm*Warm)))



# -------------------------------------------------------------------------
# OVERLAP BETWEEN FISHERIES AND SURVEY ------------------------------------
# -------------------------------------------------------------------------

# 1- Time = Season --------------------------------------------------------
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


# - Average accross all years ---------------------------------------------
D_cold_gut_m <- apply(D_cold_gut,c(1,2),mean)
D_warm_gut_m <- apply(D_warm_gut,c(1,2),mean)



# Dplyr format ------------------------------------------------------------
D_warm_gut_m <- melt(D_warm_gut_m)
colnames(D_warm_gut_m) <- c("knot","season","Warm")

D_cold_gut_m <- melt(D_cold_gut_m)
colnames(D_cold_gut_m) <- c("knot","season","Cold")

D_f <- left_join(D_warm_gut_m,D_cold_gut_m)




# - survey CPUE -----------------------------------------------------------
load("04_outputs/EBS_grid/EnvCxSeason/survey/fit.RData")
fit_s <- fit
dim(fit_s$Report$D_gct)

D_s <- melt(fit_s$Report$D_gct)
colnames(D_s) <- c("knot","cat","year", "value")
D_s <- D_s %>% mutate(Year= year+1981)

D_s_cold <- D_s %>% filter(Year %in% cold_year)
D_s_cold_m <- D_s_cold %>% group_by(knot) %>% summarise(Cold_s=mean(value))

D_s_warm <- D_s %>% filter(Year %in% warm_year)
D_s_warm_m <- D_s_warm %>% group_by(knot) %>% summarise(Warm_s=mean(value))

D_s_m <- left_join(D_s_warm_m,D_s_cold_m)


# - Join Survey and fisheries ---------------------------------------------

D <- left_join(D_f ,D_s_m)

O_season_temp <- D %>% mutate(prop_cold_s =Cold_s /(Cold_s+Cold),prop_cold_f = Cold/(Cold_s+Cold_s),
                              prop_warm_s =Warm_s /(Warm_s+Warm),prop_warm_f = Warm/(Warm_s+Warm)) %>% 
                 mutate(num_cold=prop_cold_s*prop_cold_f,num_warm= prop_warm_s*prop_warm_f)

O_season <- O_season_temp %>% group_by(season) %>% summarise(O_cold = (2*sum(num_cold))/(sum(prop_cold_s*prop_cold_s)+sum(prop_cold_f*prop_cold_f)),
                                                             O_warm = (2*sum(num_warm))/(sum(prop_warm_f*prop_warm_f)+sum(prop_warm_s*prop_warm_s)))


temp <- c(O_season$O_cold,O_season$O_warm)
temp2 <- cbind(as.numeric(temp),c("early","int","late"),c("Cold","Cold","Cold","Warm","Warm","Warm"))
temp3 <- as.data.frame(temp2)
colnames(temp3)<- c("Overlap","Season","EnvCond")
temp3 <- as_tibble(temp3)
temp3$Overlap <- as.numeric(temp3$Overlap)

p <- ggplot()+ geom_point(data=temp3 ,aes(x=Season, y=Overlap,color=EnvCond),size=3)+theme_bw()
p

png(paste(workdir,'/Overlap_averaged','.png',sep=''), height =14 , width = 25, units = 'cm', res=500)#, res=500)
p
dev.off()
p

# 2- Time = Year --------------------------------------------------------
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

O_season_temp <- D %>% mutate(num=value*value_s)

O_year <- O_season_temp %>% group_by(season,Year,EnvCond) %>% summarise(O = (2*sum(num))/(sum(value*value)+sum(value_s*value_s)))


# - Prop ------------------------------------------------------------------

D_prop <- D %>% mutate(prop_f = value/(value_s+value), prop_s= value_s/(value_s+value))


O_season_temp <-D_prop %>% mutate(num=prop_f*prop_s,den1=prop_f*prop_f,den2=prop_s*prop_s )
O_year <- O_season_temp %>% group_by(season,Year,EnvCond) %>% summarise(O = (2*sum(num))/(sum(den1)+sum(den2)))


O_year <- O_season_temp %>% group_by(season,Year,EnvCond) %>% summarise(O = (2*sum(num))/(sum(prop_f*prop_f)+sum(prop_s*prop_s)))


p <- ggplot(O_year)+
  geom_line(aes(x=Year,y=O,color=as.factor(season)))+geom_point(aes(x=Year, y=O,shape=EnvCond),size=3)+theme_bw()
p

png(paste(workdir,'/Overlap_seasons','.png',sep=''), height =14 , width = 25, units = 'cm', res=500)#, res=500)
p
dev.off()


p <- ggplot(O_year %>% filter(season=="2"))+
  geom_line(aes(x=Year,y=O,color=as.factor(season)))+geom_point(aes(x=Year, y=O,shape=EnvCond),size=3)+theme_bw()
p

png(paste(workdir,'/Overlap_intermediate','.png',sep=''), height =14 , width = 25, units = 'cm', res=500)#, res=500)
p
dev.off()






p.resol_ab <- #ggplot() +
  ggplot( data=bathy.dat3,aes(x = lon, y = lat ))+ geom_contour(aes(z = depth),breaks  = c(-100,-200,-50),col="brown") +
  # geom_point(data= fit$data_frame,aes(x=E_km, y= N_km),col="grey",size=3,shape=3)+
  geom_point(data= as_tibble(COG2_ab),aes(x= x, y=y
                                          ,
                                          col=Env.Cond#,
                                          #shape=seasons
  )
  ,size=4)+theme_bw()+
  scale_colour_manual(name = "Env.Cond", values = c("Cold"="blue","Warm"="red"))+
  geom_polygon(data=map_data("world") , aes(x=long, y = lat, group = group), fill="grey")+coord_fixed(xlim = xlims,ylim=ylims) +
  
  #stat_contour( data=bathy.dat,aes(x = lon, y = lat,z = depth), breaks  = c(-200,-100,-50),col="black",size=0.2)+
  geom_text_contour(aes(z = depth),breaks  = c(-200,0 ,-100,10,-50),, stroke = 0.2,col="brown")

p.resol_ab + facet_wrap(~seasons,ncol=1)



png(paste(workdir,'/COG_seasons','.png',sep=''), height =14 , width = 25, units = 'cm', res=500)#, res=500)
p.resol_ab + facet_wrap(~seasons,ncol=1)
dev.off()




# LIC ---------------------------------------------------------------------

D_temp <- D %>% group_by(Year,EnvCond,season) %>% summarise(LIC=sum(value*value_s)/sqrt(sum(value_s*value_s)*sum(value*value)))



p <- ggplot(D_temp)+
  geom_line(aes(x=Year,y=LIC,color=as.factor(season)))+geom_point(aes(x=Year, y=LIC,shape=EnvCond),size=3)+theme_bw()
p

png(paste(workdir,'/LIC_seasons','.png',sep=''), height =14 , width = 25, units = 'cm', res=500)#, res=500)
p
dev.off()




