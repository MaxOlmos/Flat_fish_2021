library(tidyverse)
library(ggthemes)
library(lubridate)
df <- read_csv("yfs_length.csv")
#AKFIN File had crap names so re-did
names(df) <- c("Year","T_Table","Haul_Join","Port_Join","Cruise","Permit","Length_Seq","Haul_Offload_Date","Haul_Offload","Fishing_DepthFathom","Bottom_Depth_Fathoms","IFQ_Flag","Gear","Gear_Description","Performance","Performance_Description","FMP_Area","FMP_Subarea","NMFS_Area","Species_Code","Species_Name","Sex","Length","Frequency","Eggs","Viability","Injury","Sample_System","Vessel","Vessel_Length","LatDD_Start","LonDD_Start","LatDD_End","LonDD_End","Received_from_NORPAC","Loaded_to_Repository")

glimpse(df)
# Make a set for maxime
dfm <- df %>% transmute(Year, Frequency, Length, Sex,date=dmy(Haul_Offload_Date),month=month(date),depth=Bottom_Depth_Fathoms/1.8288, quarter=quarter(date),
	                       NMFS_Area,
	                       LonDD_Start, LonDD_End,
	                       LatDD_Start, LatDD_End,
	                       ) 

write_csv(dfm,"For_Max.csv")
library(ggplot2)
dfm %>% group_by(Year,Sex) %>% filter(Sex!="U") %>% summarize(N=sum(Frequency),mnlen=sum(Length*Frequency)/N) %>% ggplot(aes(x=Year,y=mnlen,color=Sex)) + 
                                                 geom_line(size=2) + ylab("Mean Length (cm)") + theme_few()

dfs <- dfm %>% filter(Sex!="U") %>% transmute(Year = Year, Frequency, Length, Sex,date,month, quarter) %>% 
           group_by(Year,Sex,quarter) %>% summarize(N=sum(Frequency),mnlen=sum(Length*Frequency)/N) 

dfs %>% mutate(qf=as.factor(quarter)) %>% ggplot(aes(x=Year,y=N,fill=qf,shape=Sex,color=qf)) + geom_line(size=1) + geom_point(size=4) + ylab("Number of samples ") + facet_grid(quarter~.) + theme_few()
dfs %>% filter(Year>1990) %>% ggplot(aes(x=Year,y=N,shape=Sex)) + geom_line(color="blue") + geom_point(size=3,color="salmon") + ylab("Number of samples ") + facet_grid(quarter~.) + theme_few()
