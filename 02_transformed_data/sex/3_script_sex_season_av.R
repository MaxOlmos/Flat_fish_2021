
###########################################################################
########           Script CPUE VAST sex
###########################################################################

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
rm(list=ls())

# Load package ------------------------------------------------------------
# Install FishStatsUtils from CRAN
library(tidyverse)
library(ggthemes)
library(lubridate)
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
library(viridis)
#devtools::install_github("mnbram/gggibbous")
#https://cran.r-project.org/web/packages/gggibbous/vignettes/gggibbous.html
library(gggibbous)

getwd()
DateFile <- "02_transformed_data/sex/"
#setwd(DateFile)


# PREPARE DATA ----------------------------------------------------------
# -------------------------------------------------------------------------
getwd()
dfm <- read_csv(paste0(DateFile,"For_Max.csv"))
user_region <- readRDS('02_transformed_data/user_grid_outEBS/user_region.rds')
plot(user_region$Lon,user_region$Lat)

# Script Jim --------------------------------------------------------------
dfm %>% group_by(Year,Sex) %>% filter(Sex!="U",Year>1999, Year<2020) %>% summarize(N=sum(Frequency),mnlen=sum(Length*Frequency)/N) %>% ggplot(aes(x=Year,y=mnlen,color=Sex)) + 
  geom_line(size=2) + ylab("Mean Length (cm)") + theme_few()

dfs <- dfm %>% filter(Sex!="U",Year>2000, Year<2019) %>% transmute(Year = Year, Frequency, Length, Sex,date,month, quarter) %>% 
  group_by(Year,Sex,quarter) %>% summarize(N=sum(Frequency),mnlen=sum(Length*Frequency)/N) 

dfs %>% mutate(qf=as.factor(quarter)) %>% ggplot(aes(x=Year,y=N,fill=qf,shape=Sex,color=qf)) + geom_line(size=1) + geom_point(size=4) + ylab("Number of samples ") + facet_grid(quarter~.) + theme_few()
dfs %>% filter(Year>1990) %>% ggplot(aes(x=Year,y=N,shape=Sex)) + geom_line(color="blue") + geom_point(size=3,color="salmon") + ylab("Number of samples ") + facet_grid(quarter~.) + theme_few()


# Select variable of interest ---------------------------------------------
# -------------------------------------------------------------------------
dfm$LonDD_End <- dfm$LonDD_End %>% replace_na(0)
dfm$LonDD_Start <- dfm$LonDD_Start %>% replace_na(0)
dfm$LatDD_End <- dfm$LatDD_End %>% replace_na(0)
dfm$LatDD_Start <- dfm$LatDD_Start %>% replace_na(0)

df <- dfm %>% dplyr::filter(Sex!="U",Year>2000, Year<2020) %>% dplyr::mutate(lattemp=(LatDD_End+LatDD_Start), lontemp=(LonDD_End+LonDD_Start))
df <- df %>% dplyr::mutate(lat=ifelse(lattemp>90,lattemp/2,lattemp),lon=ifelse(lontemp< -181,lontemp/2,lontemp))
df %>% filter(Year==2019) #== 157 observation does not have spatial coordinates


# NA observation ----------------------------------------------------------
df %>% filter(lat==0) #== 157 observation does not have spatial coordinates


# Select only observation with spatial coordinates ------------------------
# -------------------------------------------------------------------------
df <- df %>% filter(lat>0,lon<0)
df %>% filter(lon==0)
df %>% filter(lon > -159)

# Merge sex dat wiht ADFG grid --------------------------------------------
# -------------------------------------------------------------------------

lat <- seq(53.75,62,0.5)
long <-seq(-179.5,-156.5,1)
ADFG <- expand.grid(lat, long, stringsAsFactors = FALSE)
colnames(ADFG) <-c("lat","long") 
plot(ADFG)

NN_catch  = RANN::nn2(data = ADFG[, c("long", "lat")],
                      query = df[, c("lon", "lat")],
                      k = 1)

ADFG[NN_catch$nn.idx[,1],c("long","lat")]
df2 <- cbind(df[,c("Year","Frequency", "Length", "Sex","date","month","depth")],ADFG[NN_catch$nn.idx[,1],c("long","lat")])

# Select week -------------------------------------------------------------

week <- as.character(df2$date, format="%U")
df3 <- cbind(df2,week)

df4 <- df3 %>% dplyr::group_by(Year,Sex, long, lat,week) %>% dplyr::summarise(N=sum(Frequency))
df4b <- df3 %>% dplyr::group_by(Year, long, lat,week) %>% dplyr::summarise(Ntot=sum(Frequency))

df4byear <- df3 %>% dplyr::group_by(long, lat,week) %>% dplyr::summarise(Ntotyear=sum(Frequency))

df5 <- left_join(df4,df4b)
df5_2 <-left_join(df5,df4byear)

# aaverage by sason

early <- c(12:20) # end of march- mid/end may
int <- c(21:30)#end of may to mid/end July
late<-(31:38)# 16 august to end september

Late <- cbind(rep("late",length(late)),late)
colnames(Late) <- c("season","week")
Early <- cbind(rep("early",length(early)),early)
colnames(Early) <- c("season","week")
Int <- cbind(rep("int",length(int)),int)
colnames(Int) <- c("season","week")

Season <- rbind(Early,Int,Late)

df5$week <- as.character(df5$week)
df6_temp <- right_join(as_tibble(df5), as_tibble(Season))
df6 <- df6_temp %>% dplyr::group_by(Year,season,Sex, long,lat) %>% dplyr::summarize(N_tot=sum(N), Ntot2 = sum(Ntot),Prop=N_tot/Ntot2)

df6<- df6 %>% mutate(right=ifelse(Sex=="F","TRUE","FALSE"))
df6$right <- as.logical(df6$right)



# -------------------------------------------------------------------------
# Merge fisheries data with sex data --------------------------------------
# -------------------------------------------------------------------------
load(file="02_transformed_data/observer/CPUE_catchability_adfg.RData")
CPUE_catchability <- CPUE_catchability %>% filter(year>2000)

df6F <- df6 %>% filter(Sex=='F')
df6F_0 <- df6 %>% filter(Sex=='M',Prop==1) %>% transmute(Year, Sex="F",Prop=0, long,lat,season,Ntot2=0,Ntot2)

df7F <- bind_rows(df6F ,df6F_0)

colnames(df7F) <- c( "year" ,"Season","Sex",  "long" ,"lat" , "N","Ntot","prop","right")

CPUE_catchability_s <- left_join(CPUE_catchability,df7F)

dim(CPUE_catchability)
dim(df7F)
dim(CPUE_catchability_s)


# proportion of NA --------------------------------------------------------
CPUE_catchability_s %>% filter(is.na(prop),Season=='int') # 1192 out of 4,712 observation 


# plot --------------------------------------------------------------------
# accross year
    
xlims <- xlims <- c(-177 ,-158)
ylims <- c(53, 63)
df6 <- df6 %>% mutate(Season=season,year=Year)
y <- sort(unique(df6F$Year))


  png(paste('Season_sex_CPUE1.png',sep=''), res=600,width = 25,
      height = 13,units="in")
  p <- ggplot() +
    geom_point(data=CPUE_catchability_s %>% filter(year %in% c(2001:2010)),aes(x=long,y=lat),col="green3",shape=4,size=2)  +
    geom_polygon(data=map_data("world") , aes(x=long, y = lat, group = group), fill="grey")+
    geom_moon(data=df6 %>% filter(year %in% c(2001:2010)),
              aes(x=long, y=lat, ratio = Prop,  fill = Sex, right = right),size=3.5,
              ,
              key_glyph = draw_key_full_moon
    )+
    #lims(x = c(0.5, 3.5), y = c(0.5, 3.5)) +
    # scale_size("size", range = c(5, 10), breaks = 2^(1:3)) +
    scale_fill_manual(values = c("firebrick1", "dodgerblue2")) +
    theme(legend.box = "horizontal") +
    coord_fixed(xlim = xlims,ylim=ylims) +
    theme_bw()
  q <- p + facet_grid(year~ Season)
  print(q)
  dev.off()

  png(paste('Season_sex_CPUE2.png',sep=''), res=600,width = 25,
      height = 13,units="in")
  p <- ggplot() +
    geom_point(data=CPUE_catchability_s %>% filter(year %in% c(2011:2019)),aes(x=long,y=lat),col="green3",shape=4,size=2)  +
    geom_polygon(data=map_data("world") , aes(x=long, y = lat, group = group), fill="grey")+
    geom_moon(data=df6 %>% filter(year %in% c(2011:2019)),
              aes(x=long, y=lat, ratio = Prop,  fill = Sex, right = right),size=3.5,
              ,
              key_glyph = draw_key_full_moon
    )+
    #lims(x = c(0.5, 3.5), y = c(0.5, 3.5)) +
    # scale_size("size", range = c(5, 10), breaks = 2^(1:3)) +
    scale_fill_manual(values = c("firebrick1", "dodgerblue2")) +
    theme(legend.box = "horizontal") +
    coord_fixed(xlim = xlims,ylim=ylims) +
    theme_bw()
  q <- p + facet_grid(year~ Season)
  print(q)
  dev.off()

  
  

# -------------------------------------------------------------------------
# replace NA in fihseries data with closest within the season and CP --------
# -------------------------------------------------------------------------
CPUE_new <- NULL
  CPUE_noNA %>% filter(Season=="early",Sex=="F")
for (i in c("early","int","late")){
  for (j in c("cold", "warm")){
    CPUE <- CPUE_catchability_s %>% filter(Season==i, CP==j)
    CPUE_noNA <- CPUE %>% filter(!is.na(prop))
    CPUE_NA <- CPUE %>% filter(is.na(prop))
    
    NN_temp  = RANN::nn2(data = CPUE_noNA[, c("long", "lat")],
                         query = CPUE_NA[, c("long", "lat")],
                         k = 1)
    
    CPUE_noNA[NN_temp$nn.idx[,1],c("long","lat","prop")]
    
    
    temp <- cbind(CPUE_NA[,c("year", "adfg", "CPUE.mean","CP","Season","long","lat")],
                  CPUE_noNA[NN_temp$nn.idx[,1],c("prop","Sex","N","Ntot")]
    )
    
    
    temp_new <- bind_rows(CPUE_noNA,temp)
    
    CPUE_new <- bind_rows(CPUE_new,temp_new)
  }
}
CPUE_new %>% filter(Season=="early",Sex=="F",year==2001)
df6 %>% filter(Season=="early",Sex=="F",year==2001)

# plot ------------------------------------------------------
CPUE_new2 <- CPUE_new %>% mutate(Sex="M",prop= 1-prop)
CPUE_new3 <- bind_rows(CPUE_new,CPUE_new2)
CPUE_new3<- CPUE_new3 %>% mutate(right=ifelse(Sex=="F","TRUE","FALSE"))
CPUE_new3$right <- as.logical(CPUE_new3$right)




png(paste('Season_sex_CPUE1_new.png',sep=''), res=600,width = 25,
    height = 13,units="in")
p <- ggplot() +
  geom_polygon(data=map_data("world") , aes(x=long, y = lat, group = group), fill="grey")+
  geom_moon(data=CPUE_new3 %>% filter(year %in% c(2001:2010)),
            aes(x=long, y=lat, ratio = prop,  fill = Sex, right = right),size=3.5,
            ,
            key_glyph = draw_key_full_moon
  )+
  #lims(x = c(0.5, 3.5), y = c(0.5, 3.5)) +
  # scale_size("size", range = c(5, 10), breaks = 2^(1:3)) +
  scale_fill_manual(values = c("firebrick1", "dodgerblue2")) +
  theme(legend.box = "horizontal") +
  coord_fixed(xlim = xlims,ylim=ylims) +
  theme_bw()
q <- p + facet_grid(year~ Season)
print(q)
dev.off()

png(paste('Season_sex_CPUE2_new.png',sep=''), res=600,width = 25,
    height = 13,units="in")
p <- ggplot() +
  geom_polygon(data=map_data("world") , aes(x=long, y = lat, group = group), fill="grey")+
  geom_moon(data=CPUE_new3 %>% filter(year %in% c(2011:2019)),
            aes(x=long, y=lat, ratio = prop,  fill = Sex, right = right),size=3.5,
            ,
            key_glyph = draw_key_full_moon
  )+
  #lims(x = c(0.5, 3.5), y = c(0.5, 3.5)) +
  # scale_size("size", range = c(5, 10), breaks = 2^(1:3)) +
  scale_fill_manual(values = c("firebrick1", "dodgerblue2")) +  theme(legend.box = "horizontal") +
  coord_fixed(xlim = xlims,ylim=ylims) +
  theme_bw()
q <- p + facet_grid(year~ Season)
print(q)
dev.off()


# define female/male ------------------------------------------------------

CPUE_female <- CPUE_new %>% 
  mutate(CPUE.mean = CPUE.mean*prop) %>% 
  dplyr::select(year,CPUE.mean,long,lat,CP,Season) %>% mutate(Sex="F")

CPUE_male <- CPUE_new %>% 
  mutate(CPUE.mean = CPUE.mean*(1-prop)) %>% 
  dplyr::select(year,CPUE.mean,long,lat,CP,Season) %>% mutate(Sex="M")

CPUE_sex <- bind_rows(CPUE_female,CPUE_male)
getwd()
save(CPUE_sex,file="02_transformed_data/sex/CPUE_sex.RData")

# -------------------------------------------------------------------------

