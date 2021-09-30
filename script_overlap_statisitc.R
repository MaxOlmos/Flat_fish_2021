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

# 1- Time = Season --------------------------------------------------------
# -------------------------------------------------------------------------



# - fisheries CPUE --------------------------------------------------------
load("04_outputs/EBS_grid/EnvCxSeason/fisheries/c(3,3,3,3,2,2)/fit.RData")
fit_f <- fit
dim(fit_f$Report$D_gct)

# Calculate abundance because we want proportion
Spatial_List <- fit$spatial_list
a_k <- c(Spatial_List$a_xl[, 1], rep(0,Spatial_List$MeshList$anisotropic_mesh$n-nrow(Spatial_List$MeshList$loc_x)))
knot_i <- seq(1,length(a_k),1)
a_k <- as.data.frame(cbind(a_k,knot_i))
colnames(a_k) < c("a_k","knot")

#define cold or warm 

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
      D_cold_gut[g,u,t] <- fit_f$Report$D_gct[g,1,cold[t]] * a_k[g,1] * exp(fit_f$Report$Phi1_sk[g,cold_season[u]]) 
    }
  }
}

for ( t in 1:length(warm)){
  for ( u in 1 :length(warm_season)){
    for ( g in 1:n_knot){
      D_warm_gut[g,u,t] <- fit_f$Report$D_gct[g,1,warm[t]] *a_k[g,1] * exp(fit_f$Report$Phi1_sk[g,warm_season[u]])  
    }
  }
}


# - Average accross all years ---------------------------------------------
D_cold_gut_m <- apply(D_cold_gut,c(1,2),mean)
D_warm_gut_m <- apply(D_warm_gut,c(1,2),mean)



# Dplyr format 
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
dim(D)

O_season_temp <- D %>% group_by(season) %>% summarise(prop_warm_s = Warm_s/sum(Warm_s), prop_warm_f = Warm/sum(Warm), prop_cold_s = Cold_s/sum(Cold_s), prop_cold_f = Cold/sum(Cold)) 
O_season_temp <- left_join(O_season_temp,D[,c(1,2)])


O_season <- O_season_temp %>% group_by(season) %>% summarise(O_cold = (2*sum(prop_cold_s*prop_cold_f))/(sum(prop_cold_s*prop_cold_s)+sum(prop_cold_f*prop_cold_f)),
                                                             O_warm = (2*sum(prop_warm_s*prop_warm_f))/(sum(prop_warm_f*prop_warm_f)+sum(prop_warm_s*prop_warm_s)))


temp <- c(O_season$O_cold,O_season$O_warm)
temp2 <- cbind(as.numeric(temp),c("early","int","late"),c("Cold","Cold","Cold","Warm","Warm","Warm"))
temp3 <- as.data.frame(temp2)
colnames(temp3)<- c("Overlap","Season","EnvCond")
temp3 <- as_tibble(temp3)
temp3$Overlap <- as.numeric(temp3$Overlap)

p <- ggplot(data=temp3 ,aes(x=Season, y=Overlap,color=EnvCond),size=5)+theme_bw()+geom_point(size=5)
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
      D_cold_gut[g,u,t] <- fit_f$Report$D_gct[g,1,cold[t]]  *a_k[g,1]* exp(fit_f$Report$Phi1_sk[g,cold_season[u]]) 
    }
  }
}

for ( t in 1:length(warm)){
  for ( u in 1 :length(warm_season)){
    for ( g in 1:n_knot){
      D_warm_gut[g,u,t] <- fit_f$Report$D_gct[g,1,warm[t]] *a_k[g,1] * exp(fit_f$Report$Phi1_sk[g,warm_season[u]])  
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
dim(D)

D %>% filter(year==1,season==1,EnvCond=="Cold")
O_season_temp <- D %>% group_by(Year,season,EnvCond) %>% summarise(prop_s =  value_s/sum(value_s), prop =  value/sum(value)) 
O_season_temp <- left_join(O_season_temp,D[,c(1,2)])


O_year <- O_season_temp %>% group_by(season,Year,EnvCond) %>% summarise(O = (2*sum(prop_s*prop))/(sum(prop_s*prop_s)+sum(prop*prop)))


# - Plot ------------------------------------------------------------------

O_year <- O_year %>% mutate(Seasons=recode(season, "1" ="Early","2" ="Intermediate", "3"="Late"))
p <- ggplot(O_year)+
  geom_line(aes(x=Year,y=O,color=Seasons),size=2)+geom_point(aes(x=Year, y=O,size=EnvCond))+theme_bw()
p 

png(paste(workdir,'/Overlap_seasons','.png',sep=''), height =14 , width = 25, units = 'cm', res=500)#, res=500)
p
dev.off()


# considers weighted max density ------------------------------------------
# -------------------------------------------------------------------------

# Nseason/Nmax by year and season
D <- left_join(D_f ,D_s)
dim(D)
D <- as_tibble(D)

D_temp <- D %>% group_by(Year,season) %>% summarize(N_season= sum(value)) 
D_temp2 <-D_temp %>% group_by(Year) %>% summarise(N_max=max(N_season)) 
D_temp3 <- left_join(D_temp,D_temp2)
D_temp3 %>% filter(Year==2006)
D_temp4 <- left_join(D,D_temp3) 


O_season_temp <- D_temp4 %>% group_by(Year,season,EnvCond) %>% summarise(sum_value_s=sum(value_s), sum_value=sum(value)) 
O_season_temp2 <- left_join(D_temp4,O_season_temp) %>% mutate(prop_s =  value_s/sum_value_s, prop =  value/sum_value) 


O_season_temp3 <- O_season_temp2 %>% group_by(season,Year,EnvCond) %>% summarise(num=sum(prop_s*prop),den=(sum(prop_s*prop_s)+sum(prop*prop)))
O_season_temp4 <- left_join(O_season_temp3,D_temp3)


O_year <- O_season_temp4 %>% group_by(season,Year,EnvCond) %>% summarise(O = (2*num*(N_season/N_max))/(den))

O_year <- O_year %>% mutate(Seasons=recode(season, "1" ="Early","2" ="Intermediate", "3"="Late"))
p <- ggplot(O_year)+
  geom_line(aes(x=Year,y=O,color=Seasons),size=2)+geom_point(aes(x=Year, y=O,size=EnvCond))+theme_bw()
p 



O_year <- O_season_temp2 %>% group_by(season,Year,EnvCond) %>% summarise(O = (2*(N_season/N_max)*sum(prop_s*prop))/(sum(prop_s*prop_s)+sum(prop*prop)))
O_year2 <-O_year %>% group_by(season,Year,EnvCond) %>% summarise(Om=mean(O)) 
O_year2 <- O_year2 %>% mutate(Seasons=recode(season, "1" ="Early","2" ="Intermediate", "3"="Late"))
p <- ggplot(O_year2)+
  geom_line(aes(x=Year,y=Om,color=Seasons),size=2)+geom_point(aes(x=Year, y=Om,size=EnvCond))+theme_bw()
p 


# Nseason/Nmax by knot and year season
D <- left_join(D_f ,D_s)
dim(D)
D <- as_tibble(D)

D_temp <- D %>% group_by(Year,season,knot) %>% summarize(N_season= sum(value)) 
D_temp2 <-D_temp %>% group_by(Year,knot) %>% summarise(N_max=max(N_season)) 
D_temp3 <- left_join(D_temp,D_temp2)
D_temp3 %>% filter(Year==2006)
D_temp4 <- left_join(D,D_temp3) 


O_season_temp <- D_temp4 %>% group_by(Year,season,EnvCond) %>% summarise(sum_value_s=sum(value_s), sum_value=sum(value)) 
O_season_temp2 <- left_join(D_temp4,O_season_temp) %>% mutate(prop_s =  value_s/sum_value_s, prop =  value/sum_value) 


O_season_temp3 <- O_season_temp2 %>% group_by(season,Year,EnvCond) %>% summarise(num=sum((N_season/N_max)*prop_s*prop),den=(sum(prop_s*prop_s)+sum(prop*prop)))
O_season_temp4 <- left_join(O_season_temp3,D_temp3)


O_year <- O_season_temp4 %>% group_by(season,Year,EnvCond) %>% summarise(O = (2*num)/(den))

O_year <- O_year %>% mutate(Seasons=recode(season, "1" ="Early","2" ="Intermediate", "3"="Late"))
p <- ggplot(O_year)+
  geom_line(aes(x=Year,y=O,color=Seasons),size=2)+geom_point(aes(x=Year, y=O,size=EnvCond))+theme_bw()
p 



O_year <- O_season_temp2 %>% group_by(season,Year,EnvCond) %>% summarise(O = (2*(N_season/N_max)*sum(prop_s*prop))/(sum(prop_s*prop_s)+sum(prop*prop)))
O_year2 <-O_year %>% group_by(season,Year,EnvCond) %>% summarise(Om=mean(O)) 
O_year2 <- O_year2 %>% mutate(Seasons=recode(season, "1" ="Early","2" ="Intermediate", "3"="Late"))
p <- ggplot(O_year2)+
  geom_line(aes(x=Year,y=Om,color=Seasons),size=2)+geom_point(aes(x=Year, y=Om,size=EnvCond))+theme_bw()
p 


