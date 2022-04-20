
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

df <- dfm %>% dplyr::filter(Sex!="U",Year>2000, Year<2019) %>% dplyr::mutate(lattemp=(LatDD_End+LatDD_Start), lontemp=(LonDD_End+LonDD_Start))
df <- df %>% dplyr::mutate(lat=ifelse(lattemp>90,lattemp/2,lattemp),lon=ifelse(lontemp< -181,lontemp/2,lontemp))


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
df6 <- df5 %>% mutate(prop=N/Ntot)
df6 %>% filter(Year==2003,week==21,lat<56)

df6year <- df5_2 %>% dplyr::group_by(week,Sex, long,lat) %>% dplyr::summarize(Prop =sum(N)/Ntotyear)
df6year %>% filter(Sex=="M")

df6year2 <- df6year %>% dplyr::group_by(week, Sex,long,lat) %>% dplyr::summarise(Prop=mean(Prop))

df6year3 <- df6year2 %>% mutate(right=ifelse(Sex=="F","TRUE","FALSE"))
df6<- df6 %>% mutate(right=ifelse(Sex=="F","TRUE","FALSE"))
df6$right <- as.logical(df6$right)


# plot --------------------------------------------------------------------
# -------------------------------------------------------------------------

xlims <- xlims <- c(-177 ,-158)
ylims <- c(53, 63)

# mean of prop over years -------------------------------------------------

plot <- df6year3 %>% filter(week>11, week<39)
plot %>% filter(week==37,lat==57.75,long== -165.5)
df6year3$right <- as.logical(df6year3$right)

p <- ggplot() +
  
  geom_polygon(data=map_data("world") , aes(x=long, y = lat, group = group), fill="grey")+
  geom_moon(data=df6year3 %>% filter(week>11, week<39),
    aes(x=long, y=lat, ratio = Prop, right = right, fill = Sex, size=0.01),
    key_glyph = draw_key_full_moon
  )+
  #lims(x = c(0.5, 3.5), y = c(0.5, 3.5)) +
 # scale_size("size", range = c(5, 10), breaks = 2^(1:3)) +
  scale_fill_manual(values = c("firebrick1", "dodgerblue2")) +
  theme(legend.box = "horizontal") +
  coord_fixed(xlim = xlims,ylim=ylims) +
  theme_bw()

q <- p + facet_wrap(~ week)
getwd()
png(paste('Week_sex.png',sep=''), res=600,width = 25,
    height = 13,units="in")
q
#ggsave("Week_sex.png",width = 20,
 #      height = 15,units="cm",dpi=600)

dev.off()



# -------------------------------------------------------------------------
# Merge fisheries data with sex data --------------------------------------
# -------------------------------------------------------------------------
load(file="02_transformed_data/observer/CPUE_catchability_adfg.RData")

df6F <- df6 %>% filter(Sex=='F')
df6F_0 <- df6 %>% filter(Sex=='M',prop==1) %>% transmute(Year, Sex="F",prop=0, long,lat,week,N=0,Ntot)

df7F <- bind_rows(df6F ,df6F_0)

colnames(df7F) <- c( "year" ,"Sex",  "long" ,"lat" , "week", "N","Ntot","prop","right")

CPUE_catchability_s <- left_join(CPUE_catchability,df7F)

dim(CPUE_catchability)
dim(CPUE_catchability_s)


# proportion of NA --------------------------------------------------------
CPUE_catchability_s %>% filter(is.na(prop),Season=='int') # 1192 out of 4,712 observation 


# Localisation of NA ------------------------------------------------------
p <- ggplot() +

  geom_polygon(data=map_data("world") , aes(x=long, y = lat, group = group), fill="grey")+
  geom_moon(data=df6year3 %>% filter(week>11, week<39),
            aes(x=long, y=lat, ratio = Prop, right = right, fill = Sex, size=0.01),
            key_glyph = draw_key_full_moon
  )+
  #lims(x = c(0.5, 3.5), y = c(0.5, 3.5)) +
  # scale_size("size", range = c(5, 10), breaks = 2^(1:3)) +
  scale_fill_manual(values = c("firebrick1", "dodgerblue2")) +
  theme(legend.box = "horizontal") +
  coord_fixed(xlim = xlims,ylim=ylims) +
  theme_bw()+
  geom_point(data=CPUE_catchability_s,aes(x=long,y=lat),col="green3",shape=4)

q <- p + facet_wrap(~ week)
getwd()
png(paste('Week_sex_CPUE.png',sep=''), res=600,width = 25,
    height = 13,units="in")
q

dev.off()

# accross year
y <- sort(unique(df6F$year))
plot <- df6 %>% filter(week>11, week<39)
for (i in 1:length(Year)) {
  png(paste('Week_sex_CPUE',y[i],'.png',sep=''), res=600,width = 25,
      height = 13,units="in")
p <- ggplot() +
  geom_point(data=CPUE_catchability_s %>% filter(year == y[i]),aes(x=long,y=lat),col="green3",shape=4)  +
  geom_polygon(data=map_data("world") , aes(x=long, y = lat, group = group), fill="grey")+
  geom_moon(data=plot %>% filter(Year == y[i]),
            aes(x=long, y=lat, ratio = prop,  fill = Sex, right = right, size=0.5)
            ,
            key_glyph = draw_key_full_moon
  )+
  #lims(x = c(0.5, 3.5), y = c(0.5, 3.5)) +
  # scale_size("size", range = c(5, 10), breaks = 2^(1:3)) +
  scale_fill_manual(values = c("firebrick1", "dodgerblue2")) +
  theme(legend.box = "horizontal") +
  coord_fixed(xlim = xlims,ylim=ylims) +
  theme_bw()
q <- p + facet_wrap(~ week)
print(q)
dev.off()
}


x <- df6 %>%dplyr:: filter(Year==2003,week==21,lat==54.75)
sort(x $lat)

ggplot() +
  # geom_polygon(data=map_data("world") , aes(x=long, y = lat, group = group), fill="grey")+
  geom_moon(data=plot %>% filter(Year == y[i],week==21),
            aes(x=long, y=lat, ratio = prop,  fill = Sex, right = right, size=0.5)
            ,
            key_glyph = draw_key_full_moon
  )


# -------------------------------------------------------------------------
# replace NA in fihseries data with closest within the season and CP --------
# -------------------------------------------------------------------------
CPUE_new <- NULL

for (i in c("early","int","late")){
  for (j in c("cold", "warm")){
    
CPUE <- CPUE_catchability_s %>% filter(Season==i, CP==j)
CPUE_noNA <- CPUE %>% filter(!is.na(prop))
CPUE_NA <- CPUE %>% filter(is.na(prop))

NN_temp  = RANN::nn2(data = CPUE_noNA[, c("long", "lat")],
                      query = CPUE_NA[, c("long", "lat")],
                      k = 1)

CPUE_noNA[NN_temp$nn.idx[,1],c("long","lat","prop")]


temp <- cbind(CPUE_NA[,c("year","week", "adfg", "CPUE.mean","CP","Season","long","lat")],
              CPUE_noNA[NN_temp$nn.idx[,1],c("prop","Sex","N","Ntot")]
              )


temp_new <- bind_rows(CPUE_noNA,temp)

CPUE_new <- bind_rows(CPUE_new,temp_new)
}
}
CPUE_new %>% filter(year==2001,week==12)


# plot ------------------------------------------------------
CPUE_new2 <- CPUE_new %>% mutate(Sex="M",prop=1-prop)
CPUE_new3 <- bind_rows(CPUE_new,CPUE_new2)
CPUE_new3<- CPUE_new3 %>% mutate(right=ifelse(Sex=="F","TRUE","FALSE"))
CPUE_new3$right <- as.logical(CPUE_new3$right)

p <- ggplot() +
  
  geom_polygon(data=map_data("world") , aes(x=long, y = lat, group = group), fill="grey")+
  geom_moon(data=CPUE_new3 %>% filter(week>11, week<39),
            aes(x=long, y=lat, ratio = prop, right = right, fill = Sex, size=0.01),
            key_glyph = draw_key_full_moon
  )+
  #lims(x = c(0.5, 3.5), y = c(0.5, 3.5)) +
  # scale_size("size", range = c(5, 10), breaks = 2^(1:3)) +
  scale_fill_manual(values = c("firebrick1", "dodgerblue2")) +
  theme(legend.box = "horizontal") +
  coord_fixed(xlim = xlims,ylim=ylims) +
  theme_bw()+
  geom_point(data=CPUE_catchability_s,aes(x=long,y=lat),col="green3",shape=4)

q <- p + facet_wrap(~ week)
getwd()
png(paste('Week_sex_CPUEnew2.png',sep=''), res=600,width = 25,
    height = 13,units="in")
q

dev.off()

# plot warm years
p <- ggplot() +
  
  geom_polygon(data=map_data("world") , aes(x=long, y = lat, group = group), fill="grey")+
  geom_moon(data=CPUE_new3 %>% filter(CP=="warm",week>11, week<39),
            aes(x=long, y=lat, ratio = prop, right = right, fill = Sex, size=0.01),
            key_glyph = draw_key_full_moon
  )+
  #lims(x = c(0.5, 3.5), y = c(0.5, 3.5)) +
  # scale_size("size", range = c(5, 10), breaks = 2^(1:3)) +
  scale_fill_manual(values = c("firebrick1", "dodgerblue2")) +
  theme(legend.box = "horizontal") +
  coord_fixed(xlim = xlims,ylim=ylims) +
  theme_bw()+
  geom_point(data=CPUE_catchability_s,aes(x=long,y=lat),col="green3",shape=4)

q <- p + facet_wrap(~ week)
getwd()
png(paste('Week_sex_CPUEnew_warm.png',sep=''), res=600,width = 25,
    height = 13,units="in")
q

dev.off()


#plot cold years
p <- ggplot() +
  
  geom_polygon(data=map_data("world") , aes(x=long, y = lat, group = group), fill="grey")+
  geom_moon(data=CPUE_new3 %>% filter(CP=="cold",week>11, week<39),
            aes(x=long, y=lat, ratio = prop, right = right, fill = Sex, size=0.01),
            key_glyph = draw_key_full_moon
  )+
  #lims(x = c(0.5, 3.5), y = c(0.5, 3.5)) +
  # scale_size("size", range = c(5, 10), breaks = 2^(1:3)) +
  scale_fill_manual(values = c("firebrick1", "dodgerblue2")) +
  theme(legend.box = "horizontal") +
  coord_fixed(xlim = xlims,ylim=ylims) +
  theme_bw()+
  geom_point(data=CPUE_catchability_s,aes(x=long,y=lat),col="green3",shape=4)

q <- p + facet_wrap(~ week)
getwd()
png(paste('Week_sex_CPUEnew_cold.png',sep=''), res=600,width = 25,
    height = 13,units="in")
q

dev.off()


CPUE_female <- CPUE_new %>% 
              mutate(CPUE.mean = CPUE.mean*prop) %>% 
              dplyr::select(year,week,CPUE.mean,long,lat,CP,Season) %>% mutate(Sex="F")

CPUE_male <- CPUE_new %>% 
  mutate(CPUE.mean = CPUE.mean*(1-prop)) %>% 
  dplyr::select(year,week,CPUE.mean,long,lat,CP,Season) %>% mutate(Sex="M")

CPUE_sex <- bind_rows(CPUE_female,CPUE_male)
getwd()
save(CPUE_sex,file="02_transformed_data/sex/CPUE_sex.RData")

# -------------------------------------------------------------------------

CPUE_new_F <- CPUE_new3 %>% filter(Sex=="F")

CPUE_F_cold <- CPUE_new_F %>% filter(CP=="cold") %>% group_by(week,lat,long,Season) %>% summarise(propc=sum(N)/sum(Ntot))

CPUE_F_warm <- CPUE_new_F %>% filter(CP=="warm")  %>% group_by(week,lat,long,Season) %>% summarise(propw=sum(N)/sum(Ntot))


CPUE_F <- full_join(CPUE_F_cold,CPUE_F_warm,
                    by = c( "week", "long", "lat","Season"))
#CPUE_F <- CPUE_F %>% mutate(propc = replace_na(propc,0),propw = replace_na(propw,0))

CPUE_delta <- CPUE_F %>% mutate(Dir=ifelse(is.na(propc), "warm" ,ifelse(is.na(propw),"cold" ,NA))) %>% 
              mutate(COLDvsWARM=(propc-propw)/propc) %>% 
              mutate(COLDvsWARM=ifelse(is.na(propc), propw ,ifelse(is.na(propw), propc ,COLDvsWARM)))%>% 
              mutate(Legend = case_when(is.na(Dir)&COLDvsWARM<0 ~ "<0", is.na(Dir)&COLDvsWARM==0 ~ "0",is.na(Dir)&COLDvsWARM>0 ~">0",
                                        Dir=="warm"~"WarmYears",Dir=="cold"~"ColdYears")
                     )

getwd()
bathy.dat<-read.table('01_data/bathyEBS/BeringDepth.txt',sep='') 
names(bathy.dat)<-c('lon','lat','depth')
bathy.dat$depth[bathy.dat$depth>0]<-NA#Avoid points above water

p <- ggplot() +
  
  geom_polygon(data=map_data("world") , aes(x=long, y = lat, group = group), fill="grey")+
  geom_point(data = CPUE_delta, aes(
    x = long,
    y = lat,
    color = log(sqrt(COLDvsWARM^2)+0.1),
    shape=Legend),size =2 )  +
  scale_color_viridis(name="delta(log scale)",na.value = "gray86")+
  #scale_shape_manual()+# values=c("<0"=17,">0"=1,"0"=3))+
  coord_fixed(xlim = xlims,ylim=ylims) +
  stat_contour( data=bathy.dat,aes(x = lon, y = lat,z = depth,alpha=0.5),breaks  = c(-200,-100,-50),col="brown")+
  theme_bw()

q <- p + facet_wrap(~ week)
getwd()
png(paste('COLDvsWARM.png',sep=''), res=600,width = 25,
    height = 13,units="in")
q

dev.off()

unique(CPUE_sex$week)


# Comparison --------------------------------------------------------------
# -------------------------------------------------------------------------
# male vs female
library(dplyr)

# Load fits ---------------------------------------------------------------

work_dir <- "04_outputs/user_grid_outEBS/sex/male/no2000/"
load(paste0(work_dir,"fit.RData"))
fit_male <-fit

work_dir <- "04_outputs/user_grid_outEBS/sex/female/no2000/"
load(paste0(work_dir,"fit.RData"))
fit_female <-fit

Fit <- bind_rows(fit_male,fit_female)

# index1
sd1 <- (summary(fit_female$parameter_estimates[["SD"]])[,2])
sd1 <- as.vector(sd1[which(names(sd1)=="Index_ctl")])

index1 <- as.data.frame(cbind(fit_female$Report$Index_ctl, c(2001:2019),sd1))
index1 <- index1 %>% mutate(Model="female")
colnames(index1) <- c("Index", "Years","sd","sex")
index1 <- as_tibble(index1) %>% dplyr::mutate(low=Index-0.67*sd, high=Index+0.67*sd)


# index2
sd2 <- (summary(fit_male$parameter_estimates[["SD"]])[,2])
sd2 <- as.vector(sd2[which(names(sd2)=="Index_ctl")])

index2 <- as.data.frame(cbind(fit_male$Report$Index_ctl, c(2001:2019),sd2))
index2 <- index2 %>% mutate(Model="male")
colnames(index2) <- c("Index", "Years","sd","sex")

index2 <- as_tibble(index2) %>% dplyr::mutate(low=Index-0.67*sd, high=Index+0.67*sd)

Index <- bind_rows(index1, index2)

q <- ggplot() +
  # accrross time
  geom_line(data=Index,aes(x=Years,y=Index,color=sex),size=2)+
  geom_ribbon(data=Index,aes(x=Years,ymin=low,ymax=high,fill=sex),alpha=0.2)
q

ggsave(paste0(work_dir,'Comp_male_female.png'), q)



