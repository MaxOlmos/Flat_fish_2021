###########################################################################
##### OBSERVER DATA EXPLORATION -------------------------------------------
###########################################################################
# Max_Olmos ---------------------------------------------------------------
# 8.26.2020 ---------------------------------------------------------------

rm(list=ls())

# Function to transform ADFG cells number into Lat and long ---------------
# -------------------------------------------------------------------------

#The ADFG stat 6 area codes can be broken into multiple components
#the first 2 digits can be used to determine the longitude of the ADFG box (center, right, or left edge)
#the next 3 digits can be used to determine the latitude of the ADFG box (center, top, or bottom edge)
#the last digit describes if the box has been broken into multiple sub areas(0 for no, other value for yes)
#This function takes a vector of ADFG areas and converts the areas into lat/lon coordinates that are in the center of the box

getLatLong=function(ADFGvec){
  ADFGlong=as.numeric(substring(as.character(ADFGvec), 1,2))
  ADFGlat=as.numeric(substring(as.character(ADFGvec), 3,5))
  long=-((ADFGlong+100)+0.5)
  lat=ceiling(2*ADFGlat/10)/2+0.25
  return(cbind(long,lat))
}

#an example with this function
ADFGcodes=c(695500,725430)
getLatLong(ADFGcodes)

#This function takes a vector of lats and a vector of longs (paired) and converts them into an ADFG area code
#Note that this process assigns ADFG boxes of the same size, so ADFG boxes that are broken into multiple sub areas are grouped together

getADFG=function(long,lat){
  areaLong=(-trunc(long)-100)
  areaLat1=trunc(lat)
  areaLat2=lat-trunc(lat)
  areaLat2[areaLat2<0.5]=0
  areaLat2[areaLat2>=0.5]=3
  ADFG=areaLong*10000+areaLat1*100+areaLat2*10
  return(ADFG)
}

#an example
longs=c(-169.5,-172.5)
lats=c(55.25,54.75)
getADFG(longs,lats)

# Load package ------------------------------------------------------------
# -------------------------------------------------------------------------

library(readr)
library(tidyr)
library(tidyverse)
library(dplyr)
library(tibble)
library(ggplot2)
library(sf)
library(mapproj)
library(maps)
library(viridis)
library(FishData)
library(reshape)
library(lubridate)
library(grid)
getwd()
# Load Data ---------------------------------------------------------------
# -------------------------------------------------------------------------
#Units are tons. 
#obstot is the total catch in that tow
#(flats plus all other stuff, and pflat is proportion flats (the 3 spp) in the catch.
# Effort is towDur (Tow duration) within each week/area.
Data_0 <- read_csv("01_data/observer/flats.csv", col_names = TRUE)
file_output <- "02_transformed_data/observer"

# - Check NA --------------------------------------------------------------
Data_NA <- Data_0 %>% filter(is.na(lat))
(unique(Data_NA$year))
Data_NA %>% filter(year == 2015)

# - Define CPUE -----------------------------------------------------------
Data <- Data_0 

# - Determine week --------------------------------------------------------
date <- as.Date(x = paste0(Data$year, "-", Data$week, "-1"),
                format = "%Y-%U-%u")

Data_date <- cbind(Data, date)

Data_date <- Data_date %>%
  mutate(month_day = str_sub(date, start = 6, end = 10)) %>%
  mutate(month = str_sub(date, start = 6, end = 7))


# - select only summer month : from May to August -------------------------
Data_summer <-
  Data_date %>% as_tibble %>% filter(month %in% c("03","04", "05", "06", "07", "08","09")) %>% drop_na()

p <- ggplot(data=Data_summer, aes(x=week, y=log(CPUE)))+geom_point()+ 
  stat_smooth( method = "gam", colour="black",size=2, formula = y ~ s(x))
p

# -Plot years --------------------------------------------------------
p <- ggplot(data=Data_summer, aes(x=year, y=log(CPUE)))+geom_point()+ 
  stat_smooth( method = "gam", colour="black",size=2, formula = y ~ s(x))
p

png(paste0(file_output, '/logCPUE=f(years_summer)', '.png'), height =7 , width = 5, units = 'in', res=600)
p
dev.off()


# Define covariates -------------------------------------------------------
# -------------------------------------------------------------------------

# - Covariate: CP cold/warm years -----------------------------------------

cold <- cbind(c(2000,2006,2007,2008,2009,2010,2011,2012,2013,2017),rep('cold', length(c(2000,2006,2007,2008,2009,2010,2011,2012,2013,2017))))
warm <- cbind(c(2001,2002,2003,2004,2005,2014,2015,2016,2018,2019), rep('warm', length(c(2001,2002,2003,2004,2005,2014,2015,2016,2018,2019))))
#temp <- cbind(c(2000,2001,2017,2019,2020),rep("temp",5))
Temp <- as_tibble(rbind(cold,warm))
colnames(Temp)<- c('year','CP')
Temp$year <- as.numeric(Temp$year)

Data_adfg_1 <- left_join(Data_summer,Temp)
Data_adfg_1 <- Data_adfg_1 %>% filter(year>1999, year<2020)

# - Covariate Season = Early, Int , Lat -----------------------------------

#early <- c(13:19) # end of march- mid/end may
#int <- c(20:25)#end of may to mid/end July
#late<-(26:33)# 16 august to end september

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

Data_adfg_1$week <- as.character(Data_adfg_1$week)
CPUE_catchability <- right_join(as_tibble(Data_adfg_1), as_tibble(Season))

# - mean CPUE (kg)/ month/adfg --------------------------------------------------
CPUE_catchability_temp <- CPUE_catchability %>% group_by(year, season, adfg ,CP) %>% summarise(CPUE.mean = sum(yfswt*1000)/sum(towDur))

# - Calculate lon lat ADFG --------------------------------------------------
vect_lat_long <- getLatLong(CPUE_catchability_temp$adfg)

colnames(vect_lat_long) <- c("long","lat")
vect_lat_long<- as_tibble(vect_lat_long)
CPUE_catchability <- bind_cols(CPUE_catchability_temp,vect_lat_long)

# Save ---------------------------------------------------------------------
# -------------------------------------------------------------------------
getwd()
colnames(CPUE_catchability) <- c("year", "Season" , "adfg"  ,  "CP","CPUE.mean", "long","lat"  )
save(CPUE_catchability,file=paste0(file_output,"/CPUE_catchability_adfg_sum.RData"))

# Checks
CPUE_catchability %>% filter(is.na(week))
sort(unique(CPUE_catchability$week))

