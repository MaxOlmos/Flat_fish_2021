
# -------------------------------------------------------------------------
# Map hypothesis ----------------------------------------------------------
# -------------------------------------------------------------------------
library("rnaturalearth")
library("rnaturalearthdata")
library(devtools)
library(patchwork)


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


bathy.dat4 <- bathy.dat4 %>% filter(X< -150)
crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

data_fs_sf <- st_as_sf(bathy.dat4  , coords = c("X", "Y"),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
data_fs_sf_in <- st_as_sf(bathy.dat4 %>% filter(Y > 55), coords = c("X", "Y"),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

bathy_outer <- data_fs_sf %>% filter(Shelf=="Outer") %>% 
  concaveman::concaveman(concavity = 5/4) %>% mutate(Shelf="Outer")
bathy_middle <- data_fs_sf %>% filter(Shelf=="Middle") %>% 
  concaveman::concaveman(concavity = 5/4) %>% mutate(Shelf="Middle")
bathy_inner <- data_fs_sf_in %>% filter(Shelf=="Inner") %>% 
  concaveman::concaveman(concavity = 5/4)%>% mutate(Shelf="Inner")

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

ggsave(file= "02_transformed_data/map_migration/map_hyp.png")

