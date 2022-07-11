library(metR)
library(ggplot2)
library(VAST)                           # 3.8.0
library(ggplot2)                        # 2.10.0
library(dplyr)
library(tidyr)
library(concaveman)
library(sf)
rm(list=ls())
theme_set(theme_bw())

# Define map AK -----------------------------------------------------------

ak_map <- subset(map_data("world"), region=='USA' & subregion=='Alaska')

# Load Bathymetry ---------------------------------------------------------

bathy.dat<-read.table('C:/Users/Maxime/Documents/Git/Flat_fish_2021/01_data/bathyEBS/BeringDepth.txt',sep='') 
names(bathy.dat)<-c('lon','lat','depth')
head(bathy.dat)

load("02_transformed_data/bathyEBS/LL3.RData")
bathy.dat <- bathy.dat 
polygon_bathy <- st_as_sf(as.data.frame(LL3),coords = c("V1", "V2"))%>% 
  concaveman::concaveman(concavity = 5/4)
point <- st_as_sf(bathy.dat,coords = c("lon","lat"))

bathy.dat2 <- st_intersection(point,polygon_bathy$polygons )  
bathy.dat3 <- cbind(st_coordinates(bathy.dat2),bathy.dat2$depth)
bathy.dat3 <- as.data.frame(bathy.dat3)
colnames(bathy.dat3) <- c("X","Y","depth")


# Plot --------------------------------------------------------------------

p <-ggplot( data=bathy.dat3,aes(x = X, y = Y ))+ geom_contour(aes(z = depth),breaks  = c(-100,-200,-50),col="brown",size=1.2) +
  theme_bw()+
  geom_text_contour(aes(z = depth),breaks  = c(-200,0 ,-100,10,-50), stroke = 0.2,col="brown",size=6)+

  geom_polygon(data=ak_map, aes(x=long, y = lat, group = group), fill="black")+ coord_cartesian(xlim=c(-180,-155), ylim=c(54,62))

 
p

# Define areas -------------------------------------------------------------
plot(bathy.dat3$X,bathy.dat3$Y)
library(ggmap)

y.data<- c(CPUE_catchability$lat,62,54)
x.data<- c(CPUE_catchability$long,-168,-175)
plot(x.data, y.data)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(200,100,50,25),labcex=0.4,col='black',add=T)
polygon(wint1,col="blue")
polygon(wint2,col="blue")
polygon(wint3,col="blue")
polygon(feed1,col="blue")
polygon(feed2,col="blue")
polygon(sp1,col="blue")
polygon(sp2,col="blue")
polygon(sp3,col="blue")
points(spring_start,col="red",lwd=4)

wint1 <- locator()
wint1 <- cbind(wint1$x, wint1$y)
wint1 <- as.data.frame(wint1)
save(wint1,file="02_transformed_data/map_migration/wint1.RData")

load("02_transformed_data/map_migration/wint1.RData")

wint2 <- locator()
wint2 <- cbind(wint2$x, wint2$y)
wint2 <- as.data.frame(wint2)
save(wint2,file="02_transformed_data/map_migration/wint2.RData")

wint3 <- locator()
wint3 <- cbind(wint3$x, wint3$y)
wint3 <- as.data.frame(wint3)
save(wint3,file="02_transformed_data/map_migration/wint2.RData")


feed1 <- locator()
feed1 <- cbind(feed1$x, feed1$y)
feed1 <- as.data.frame(feed1)
save(feed1,file="02_transformed_data/map_migration/feed1.RData")

feed2 <- locator()
feed2 <- cbind(feed2$x, feed2$y)
feed2 <- as.data.frame(feed2)
save(feed2,file="02_transformed_data/map_migration/feed2.RData")

sp1 <- locator()
sp1 <- cbind(sp1$x, sp1$y)
sp1 <- as.data.frame(sp1)
save(sp1,file="02_transformed_data/map_migration/sp1.RData")

sp2 <- locator()
sp2 <- cbind(sp2$x, sp2$y)
sp2 <- as.data.frame(sp2)
save(sp2,file="02_transformed_data/map_migration/sp2.RData")

sp3 <- locator()
sp3 <- cbind(sp3$x, sp3$y)
sp3 <- as.data.frame(sp3)
save(sp3,file="02_transformed_data/map_migration/sp3.RData")


y.data<- c(CPUE_catchability$lat,62,54)
x.data<- c(CPUE_catchability$long,-168,-175)
plot(x.data, y.data)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(200,100,50,25),labcex=0.6,col='black',add=T)
polygon(wint1,col="blue")
polygon(wint2,col="blue")
polygon(wint3,col="blue")
polygon(feed1,col="blue")
polygon(feed2,col="blue")
polygon(sp1,col="blue")
polygon(sp2,col="blue")
polygon(sp3,col="blue")
points(spring_start,col="green",lwd=4)
points(un_start,col="brown",lwd=4)



spring_start <- locator()
save(spring_start,file="02_transformed_data/map_migration/spring_start.RData")
spring_end <- locator()
save(spring_end,file="02_transformed_data/map_migration/spring_end.RData")

summer_start <- locator()
save(summer_start,file="02_transformed_data/map_migration/summer_start.RData")
summer_end <- locator()
save(summer_end,file="02_transformed_data/map_migration/summer_end.RData")

un_start <- locator()
save(un_start,file="02_transformed_data/map_migration/un_start.RData")
un_end <- locator()
save(un_end,file="02_transformed_data/map_migration/un_end.RData")



vx_spring <- spring_end$x-spring_start$x
vy_spring <- spring_end$y-spring_start$y
d_spring=data.frame(x=spring_start$x, y=spring_start$y, vx= vx, vy= vy)

vx_summer <- summer_end$x-summer_start$x
vy_summer <- summer_end$y-summer_start$y
d_summer=data.frame(x=summer_start$x, y=summer_start$y, vx= vx_summer, vy= vy_summer)

vx_un <- un_end$x-un_start$x
vy_un <- un_end$y-un_start$y
d_un=data.frame(x=un_start$x, y= un_start$y, vx= vx_un, vy= vy_un)



# Bind migrations ---------------------------------------------------------


# - arrows ----------------------------------------------------------------
d_spring <- d_spring %>% mutate(Migration_Route = "Spring" )
d_summer <- d_summer %>% mutate(Migration_Route = "Summer" )
d_un <- d_un %>% mutate(Migration_Route = "Uncertain" )

d <- bind_rows(d_spring, d_summer, d_un)


# - polygons --------------------------------------------------------------

wint2s<-  smooth_ksmooth(as.matrix(wint2),smoothness = 1)
wint1s<-  smooth_ksmooth(as.matrix(wint1),smoothness = 1)
wint3s<-  smooth_ksmooth(as.matrix(wint3),smoothness = 1)
feed1s<-  smooth_ksmooth(as.matrix(feed1),smoothness = 1)
feed2s<-  smooth_ksmooth(as.matrix(feed2),smoothness = 1)
sp1s <-  smooth_ksmooth(as.matrix(sp1),smoothness = 1)
sp2s <-  smooth_ksmooth(as.matrix(sp2),smoothness = 1)
sp3s <-  smooth_ksmooth(as.matrix(sp3),smoothness = 1)

wint2s <- as_tibble(wint2s) %>% mutate(Areas = "Wintering")
wint1s <- as_tibble(wint1s) %>% mutate(Areas = "Wintering")
wint3s <- as_tibble(wint3s) %>% mutate(Areas = "Wintering")

feed1s <- as_tibble(feed1s) %>% mutate(Areas = "Feeding")
feed2s <- as_tibble(feed2s) %>% mutate(Areas = "Feeding")

sp1s <- as_tibble(sp1s) %>% mutate(Areas = "Spawning")
sp2s <- as_tibble(sp2s) %>% mutate(Areas = "Spawning")
sp3s <- as_tibble(sp3s) %>% mutate(Areas = "Spawning")

Areas1 <- bind_rows(wint1s,feed1s,sp1s)
Areas2 <- bind_rows(wint2s,feed2s,sp2s)
Areas3 <- bind_rows(wint3s,sp3s)

Areas1 <- bind_rows(wint2s,wint1s,wint3s,feed1s,feed2s,sp1s,sp2s,sp3s)

bathy.dat4 <- bathy.dat3 %>% mutate(Shelf = case_when ((depth >= -200) & (depth < -100) ~ "Outer"
                                                      ,(depth >= -100) & (depth < -50) ~ "Middle"
                                                      ,(depth >= - 50) & (depth < -1) ~ "Inner",TRUE ~ 'F')) %>% 
              filter(Shelf %in% c("Outer","Middle","Inner"))

bathy.dat4 <- as_tibble(bathy.dat4)
alphap <- 0.4

p + geom_polygon(data=as.data.frame(Areas1), aes(x=V1,y=V2,fill=Areas),alpha=alphap)+
  geom_polygon(data=as.data.frame(Areas2), aes(x=V1,y=V2,fill=Areas),alpha=alphap)+
  geom_polygon(data=as.data.frame(Areas3), aes(x=V1,y=V2,fill=Areas),alpha=alphap)+
  geom_segment(data=d, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy,colour=Migration_Route), arrow=arrow(type = 'closed'), size=1.3,alpha=1)+  
  scale_color_manual(breaks = c("Spring", "Summer", "Uncertain"),
                     values=c("purple", "Orange","grey" )) +
  geom_text(label="OuterShelf", x=-178.5,y=62,color = "brown",size=6)+
  geom_text(label="MiddleShelf", x=-174.5,y=62,color = "brown",size=6)+
  geom_text(label="InnerShelf", x=-167.5,y=62,color = "brown",size=6) 

ggsave(file= "02_transformed_data/map_migration/map.png")



# -------------------------------------------------------------------------
# Map hypothesis ----------------------------------------------------------
# -------------------------------------------------------------------------


# Define map AK -----------------------------------------------------------

ak_map <- subset(map_data("world"), region=='USA' & subregion=='Alaska')

# Load Bathymetry ---------------------------------------------------------

bathy.dat<-read.table('C:/Users/Maxime/Documents/Git/Flat_fish_2021/01_data/bathyEBS/BeringDepth.txt',sep='') 
names(bathy.dat)<-c('lon','lat','depth')
head(bathy.dat)

load("02_transformed_data/bathyEBS/LL3.RData")
bathy.dat <- bathy.dat 
polygon_bathy <- st_as_sf(as.data.frame(LL3),coords = c("V1", "V2"))%>% 
  concaveman::concaveman(concavity = 5/4)
point <- st_as_sf(bathy.dat,coords = c("lon","lat"))

bathy.dat2 <- st_intersection(point,polygon_bathy$polygons )  
bathy.dat3 <- cbind(st_coordinates(bathy.dat2),bathy.dat2$depth)
bathy.dat3 <- as.data.frame(bathy.dat3)
colnames(bathy.dat3) <- c("X","Y","depth")


bathy.dat4 <- bathy.dat3 %>% mutate(Shelf = case_when ((depth >= -200) & (depth < -100) ~ "Outer"
                                                       ,(depth >= -100) & (depth < -50) ~ "Middle"
                                                       ,(depth >= - 50) & (depth < -1) ~ "Inner",TRUE ~ 'F')) %>% 
  filter(Shelf %in% c("Outer","Middle","Inner"))

bathy.dat4 <- as_tibble(bathy.dat4)
alphap <- 0.4

library("rnaturalearth")
library("rnaturalearthdata")

bathy.dat4 <- bathy.dat4 %>% filter(X< -150)
crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

data_fs_sf <- st_as_sf(bathy.dat4  , coords = c("X", "Y"),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
data_fs_sf_in <- st_as_sf(bathy.dat4 %>% filter(Y > 55), coords = c("X", "Y"),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

bathy_outer <- data_fs_sf %>% filter(Shelf=="Outer") %>% 
  concaveman::concaveman(concavity = 5/4) 
bathy_middle <- data_fs_sf %>% filter(Shelf=="Middle") %>% 
  concaveman::concaveman(concavity = 5/4) 
bathy_inner <- data_fs_sf_in %>% filter(Shelf=="Inner") %>% 
  concaveman::concaveman(concavity = 5/4)

bathy_inner2 <- st_difference(bathy_inner,bathy_middle)
bathy_inner3 <- st_difference(bathy_inner2,bathy_outer) %>% mutate(Shelf="Inner")


world <- ne_countries(scale = "medium", returnclass = "sf")


bathy_sheld <- bind_rows(bathy_outer,bathy_middle,bathy_inner3)

col_h <- "lightsalmon"
col_m <- "lightgreen"
col_l<- "lightblue"

x <- c(-178,-173.5, -169)
y <- c(62.5,62.5,62.5)

Label_shelf <- cbind(bathy_sheld,x,y)

library(devtools)
#install_github("thomasp85/patchwork")
library(patchwork)

dCE <- data.frame(x=c(-173.2),y=c(58),vx=c(2),vy=c(0))
bathy_shelf_CE <- bind_rows(bathy_outer %>% mutate(Effect="High"),bathy_middle%>% mutate(Effect="Medium"),bathy_inner2%>% mutate(Effect="Low"))

qCE <-  ggplot()+
  geom_text(data=bathy_sheld,aes(x,y,label=Shelf),color = "brown")+
  geom_sf(data = bathy_shelf_CE,aes(fill= Effect),alpha=0.5,show.legend = FALSE )+
  scale_fill_manual(values=c(col_h, col_m,col_l), 
                     breaks=c("High", "Medium","Low"))+
  geom_sf(data = world,fill="black") + 
  coord_sf(xlim=c(-180,-155), ylim=c(54,63.5), expand = FALSE)+
  geom_segment(data=dCE, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy), arrow=arrow(type = 'closed'), size=3,alpha=1,col="brown")+
  ggtitle("Cold x Early")
 

bathy_shelf_CI <- bind_rows(bathy_outer %>% mutate(Effect="Medium"),bathy_middle%>% mutate(Effect="High"),bathy_inner2%>% mutate(Effect="High"))
dCI <- data.frame(x=c(-168.5,-163.5),y=c(58,58.5),vx=c(5,-5),vy=c(0,0))
qCI <-  ggplot()+
  geom_text(data=bathy_sheld,aes(x,y,label=Shelf),color = "brown")+
  geom_sf(data = bathy_shelf_CI,aes(fill= Effect),alpha=0.5,show.legend = FALSE )+
  scale_fill_manual(values=c(col_h, col_m,col_l)
                  )+
  geom_sf(data = world,fill="black") + 
  coord_sf(xlim=c(-180,-155), ylim=c(54,63.5), expand = FALSE)+
  geom_segment(data=dCI, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy), arrow=arrow(type = 'closed'), size=c(3,1),alpha=1,col="brown")+
  ggtitle("Cold x Int")

bathy_shelf_CL <- bind_rows(bathy_outer %>% mutate(Effect="Low"),bathy_middle%>% mutate(Effect="High"),bathy_inner2%>% mutate(Effect="Medium"))
dCL <- data.frame(x=c(-163.5),y=c(58.5),vx=c(-5),vy=c(0))
qCL <-  ggplot()+
  geom_text(data=bathy_sheld,aes(x,y,label=Shelf),color = "brown")+
  geom_sf(data = bathy_shelf_CL,aes(fill= Effect),alpha=0.5,show.legend = FALSE)+
  scale_fill_manual(values=c(col_h, col_m,col_l), 
                    breaks=c("High", "Medium","Low"))+
  geom_sf(data = world,fill="black") + 
  coord_sf(xlim=c(-180,-155), ylim=c(54,63.5), expand = FALSE)+
  geom_segment(data=dCL, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy), arrow=arrow(type = 'closed'), size=c(3),alpha=1,col="brown")+
  ggtitle("Cold x Late")



bathy_shelf_WE <- bind_rows(bathy_outer %>% mutate(Effect="High"),bathy_middle%>% mutate(Effect="High"),bathy_inner2%>% mutate(Effect="Low"))
dWE <- data.frame(x=c(-173.2),y=c(58),vx=c(4),vy=c(0))
qWE <-  ggplot()+
  geom_text(data=bathy_sheld,aes(x,y,label=Shelf),color = "brown")+
  geom_sf(data = bathy_shelf_WE,aes(fill= Effect),alpha=0.5,show.legend = FALSE)+
  scale_fill_manual(values=c(col_h, col_m,col_l), 
                    breaks=c("High", "Medium","Low"))+
  geom_sf(data = world,fill="black") + 
  coord_sf(xlim=c(-180,-155), ylim=c(54,63.5), expand = FALSE)+
  geom_segment(data=dWE, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy), arrow=arrow(type = 'closed'), size=c(3),alpha=1,col="brown")+
  ggtitle("Warm X Early")

bathy_shelf_WI <- bind_rows(bathy_outer %>% mutate(Effect="Low"),bathy_middle%>% mutate(Effect="High"),bathy_inner2%>% mutate(Effect="Medium"))
dWI <- data.frame(x=c(-163.5),y=c(58.5),vx=c(-5),vy=c(0))
qWI <-  ggplot()+
  geom_text(data=bathy_sheld,aes(x,y,label=Shelf),color = "brown")+
  geom_sf(data = bathy_shelf_WI,aes(fill= Effect),alpha=0.5,show.legend = FALSE)+
  scale_fill_manual(values=c(col_h, col_m,col_l), 
                    breaks=c("High", "Medium","Low"))+
  geom_sf(data = world,fill="black") + 
  coord_sf(xlim=c(-180,-155), ylim=c(54,63.5), expand = FALSE)+
  geom_segment(data=dWI, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy), arrow=arrow(type = 'closed'), size=c(3),alpha=1,col="brown")+
  ggtitle("Warm x Int")


bathy_shelf_WL <- bind_rows(bathy_outer %>% mutate(Effect="High"),bathy_middle%>% mutate(Effect="Medium"),bathy_inner2%>% mutate(Effect="Low"))
dWL <- data.frame(x=c(-169.5),y=c(58.5),vx=c(-4),vy=c(0))
qWL <-  ggplot()+
  geom_text(data=bathy_sheld,aes(x,y,label=Shelf),color = "brown")+
  geom_sf(data = bathy_shelf_WL,aes(fill= Effect),alpha=0.5)+
  scale_fill_manual(values=c(col_h, col_m,col_l), 
                    breaks=c("High", "Medium","Low"))+
  geom_sf(data = world,fill="black") + 
  coord_sf(xlim=c(-180,-155), ylim=c(54,63.5), expand = FALSE)+
  geom_segment(data=dWL, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy), arrow=arrow(type = 'closed'), size=c(3),alpha=1,col="brown")+
  ggtitle("Warm x Late")

 
(qCE | qCI | qCL)/ (qWE | qWI | qWL)

ggsave(file= "02_transformed_data/map_migration/map_hyp2.png")

