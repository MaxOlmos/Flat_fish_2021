###########################################################################
# CALCULATE COG -----------------------------------------------------------
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

# - Calculate COG for warm vs cold in each of the 3 seasons ---------------
data <- NULL
data <- D_warm_gut
unique(D_warm_gut$season)

COG_ab <- function(data,year_env){
  COG.y <- NULL
  
  temp1 <- NULL ; temp2 <- NULL;temp1.y <- NULL; temp2.y <- NULL;COG.y_temp=NULL
  
  temp1 <- data %>% dplyr::mutate(P =Ab*E_km) %>%
    dplyr::group_by(season) %>%
    dplyr::summarize(value=sum(P)) %>% 
    dplyr::select(value,season)
  temp1.y <- temp1$value
  
  temp2 <- data %>%
    dplyr::group_by(season) %>%
    dplyr::summarize(value=sum(Ab)) %>% 
    dplyr::select(value,season)
  
  temp2.y <- temp2$value
  
  #COG.y[t] <- temp1.y[t] / temp2.y[t]
  COG.y_temp <- temp1.y/temp2.y
  COG.y <- cbind(COG.y_temp,season )
  
  
  temp1 <- NULL; temp2 <- NULL; temp1.x <- NULL; temp2.x <- NULL; COG.x <- NULL
  #data$value <- as.numeric(data$value)
  temp1 <- data %>% dplyr::mutate(P = Ab*N_km ) %>% 
    dplyr::group_by(season) %>%
    dplyr::summarize(value=sum(P)) %>% 
    dplyr::select(value)
  temp1.x <- temp1$value
  
  temp2 <- data %>%
    dplyr::group_by(season) %>%
    dplyr::summarize(value=sum(Ab)) %>% 
    dplyr::select(value,season)
  temp2.x <- temp2$value
  
  
  COG.x_temp <- temp1.x/temp2.x
  COG.x <- cbind(COG.x_temp,season)
  
  COG <- left_join(as.data.frame(COG.y),as.data.frame(COG.x))
  colnames(COG)<- c("Y","seasons","X")
  return(COG)
  
}

COG_cold_ab <- COG_ab(D_cold_gut)
COG_cold_ab <- COG_cold_ab %>% mutate(Env.Cond = "Cold")

COG_warm_ab <- COG_ab(D_warm_gut)
COG_warm_ab <- COG_warm_ab %>% mutate(Env.Cond = "Warm")

COG_ab <- bind_rows(COG_cold_ab,COG_warm_ab)
COG_ab$Y <- as.numeric(COG_ab$Y)
COG_ab$X <- as.numeric(COG_ab$X)


# Plot --------------------------------------------------------------------
# add bathy

bathy.dat<-read.table('C:/Users/Maxime/Documents/Git/Flat_fish_2021/01_data/bathyEBS/BeringDepth.txt',sep='') 
names(bathy.dat)<-c('lon','lat','depth')

#bathy.dat$depth[bathy.dat$depth>0]<-NA#Avoid points above water
#head(bathy.dat)
#bathy.mat<- matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]



range(fit$data_frame$Lat_i)
range(fit$data_frame$Lon_i)
bathy.dat <- bathy.dat %>% filter(lon < -160, lon > -179,lat<63, lat>54)


utm <- COG_ab
utm1 <- data.frame(x=utm$Y,y=utm$X) 
coordinates(utm1) <- ~x+y 
class(utm1)
proj4string(utm1) <- CRS(paste0("+proj=utm +datum=WGS84 +units=km +zone=2"))#CRS("+proj=utm +datum=WGS84 +units=km +ellps=WGS84") 

utm2 <- spTransform(utm1,CRS("+proj=longlat +datum=WGS84 +no_defs"))

COG2_ab <- cbind(utm2@coords,COG_ab)



# Plot --------------------------------------------------------------------
xlims <- xlims <- c(-179 ,-158)
ylims <- c(52, 61)

load("02_transformed_data/bathyEBS/LL2.RData")
bathy.dat <- bathy.dat %>% filter(lon < -155, lon > -179,lat<61, lat>54)
polygon_bathy <- st_as_sf(as.data.frame(LL2),coords = c("V1", "V2"))%>% 
  concaveman::concaveman(concavity = 5/4)
point <- st_as_sf(bathy.dat,coords = c("lon","lat"))

bathy.dat2 <- st_intersection(point,polygon_bathy$polygons )  
bathy.dat3 <- cbind(st_coordinates(bathy.dat2),bathy.dat2$depth)
bathy.dat3 <- as.data.frame(bathy.dat3)
colnames(bathy.dat3) <- c("lon","lat","depth")



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



# Plot warm and cold spatial map -----------------------------------------
# -------------------------------------------------------------------------


MapDetails_List <- FishStatsUtils:: make_map_info( "Region"=settings$Region,spatial_list = fit$spatial_list, "NN_Extrap"=fit$spatial_list$PolygonList$NN_Extrap,                                                     
                                                   "Extrapolation_List"=fit$extrapolation_list ,fine_scale = FALSE)

shape_path <- "C:/Users/Maxime/Documents/Git/Flat_fish/Pheno-flatfish/ADFGv2/shapefile/"
coast_shapefile <- paste(shape_path, "ne_50m_land.shp", sep="")
ocean <- readOGR(coast_shapefile)

bathy.dat<-read.table('C:/Users/Maxime/Documents/Git/Flat_fish_2021/01_data/bathyEBS/BeringDepth.txt',sep='') 
names(bathy.dat)<-c('lon','lat','depth')

bathy.dat$depth[bathy.dat$depth>0]<-NA#Avoid points above water
head(bathy.dat)
bathy.mat<- matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]


dim(D_cold_gu_plot)
dim(D_warm_gut_plot)


# cold --------------------------------------------------------------------



  x <- as.data.frame(MapDetails_List$PlotDF)
  PlotDF1 = subset(x,Include==TRUE)
  
  unique(x$Include)
  PlotDF1 <- as.data.frame(PlotDF1)
  coords2 = cbind(PlotDF1$Lon,PlotDF1$Lat)
  coords2 = coords2[coords2[,1]<0,]
  P2 = SpatialPoints(coords2)
  P3 = SpatialPoints(coords2)
  
  rast <- raster(ncol=50,nrow=50)
  extent(rast) <- extent(coords2)
  col=colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                         "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))#colorRampPalette( c('blue','red', 'yellow'), bias=1)
  
  breakpoints  <-  seq(min( log(D_cold_gu_plot[PlotDF1$x2i,,])),round(max(log(D_cold_gu_plot[PlotDF1$x2i,,])),digits=1),(round(max(log(D_cold_gu_plot[PlotDF1$x2i,,])),digits=1)-min(log(D_cold_gu_plot[PlotDF1$x2i,,])))/70)
  range(log(D_cold_gu_plot[PlotDF1$x2i,,]))
  

  png(paste(workdir,"/Density_env1",'.png',sep=''), height = 8, width = 12, units = 'in', res=600)
  par(mfrow=c(9,3))
  par(mar=c(0.5, 6, 0.5, 6))
  for ( t in 1:9){
    for (k in 1:3){
    
  P2$data = log(D_cold_gu_plot[,k,t])[PlotDF1$x2i]
  rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
  
  plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(54,62),main=cold_year[t])
  
  image(rast.temp,col=col(length(breakpoints)-1),axes=TRUE,breaks=breakpoints,
        add=T,xlim=c(-180,-158),ylim=c(54,62))
  contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(200,100,50),labcex=0.4,col='black',add=T)
  #points(data_obs$long,data_obs$lat, pch=3)
  
  box()
  legend_x=c(0.1,0.2)
  legend_y=c(0.05,0.3)
  cex.legend=0.4
  xl = (1-legend_x[1])*par('usr')[1] + (legend_x[1])*par('usr')[2]
  xr = (1-legend_x[2])*par('usr')[1] + (legend_x[2])*par('usr')[2]
  yb = (1-legend_y[1])*par('usr')[3] + (legend_y[1])*par('usr')[4]
  yt = (1-legend_y[2])*par('usr')[3] + (legend_y[2])*par('usr')[4]
  if( diff(legend_y) > diff(legend_x) ){
    align = c("lt","rb")[2]
    gradient = c("x","y")[2]
  }else{
    align = c("lt","rb")[1]
    gradient = c("x","y")[1]
  }  
  
  plotrix::color.legend(xl=xl, yb=yb, xr=xr, yt=yt, legend=round(seq(min(breakpoints),max(breakpoints),length=4),2),
                        rect.col=col(length(breakpoints)-1), cex=0.6, align=align, gradient=gradient)
  
  
  
}
}
dev.off()

# warm  --------------------------------------------------------------------


x <- as.data.frame(MapDetails_List$PlotDF)
PlotDF1 = subset(x,Include==TRUE)

unique(x$Include)
PlotDF1 <- as.data.frame(PlotDF1)
coords2 = cbind(PlotDF1$Lon,PlotDF1$Lat)
coords2 = coords2[coords2[,1]<0,]
P2 = SpatialPoints(coords2)
P3 = SpatialPoints(coords2)

rast <- raster(ncol=50,nrow=50)
extent(rast) <- extent(coords2)
col=colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                       "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))#colorRampPalette( c('blue','red', 'yellow'), bias=1)

breakpoints  <-  seq(min( log(D_cold_gu_plot[PlotDF1$x2i,,])),round(max(log(D_cold_gu_plot[PlotDF1$x2i,,])),digits=1),(round(max(log(D_cold_gu_plot[PlotDF1$x2i,,])),digits=1)-min(log(D_cold_gu_plot[PlotDF1$x2i,,])))/70)
range(log(D_cold_gu_plot[PlotDF1$x2i,,]))

png(paste(workdir,"/Density_warm",'.png',sep=''), height = 8, width = 12, units = 'in', res=600)
par(mfrow=c(9,3))
par(mar=c(0.5, 6, 0.5, 6))
for ( t in 1:9){
  for (k in 1:3){
    
    P2$data = log(D_warm_gut_plot[,k,t])[PlotDF1$x2i]
    rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
    
    plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(54,62),main=warm_year[t])
    
    image(rast.temp,col=col(length(breakpoints)-1),axes=TRUE,breaks=breakpoints,
          add=T,xlim=c(-180,-158),ylim=c(54,62))
    contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(200,100,50),labcex=0.4,col='black',add=T)
    #points(data_obs$long,data_obs$lat, pch=3)
    
    box()
    legend_x=c(0.1,0.2)
    legend_y=c(0.05,0.3)
    cex.legend=0.4
    xl = (1-legend_x[1])*par('usr')[1] + (legend_x[1])*par('usr')[2]
    xr = (1-legend_x[2])*par('usr')[1] + (legend_x[2])*par('usr')[2]
    yb = (1-legend_y[1])*par('usr')[3] + (legend_y[1])*par('usr')[4]
    yt = (1-legend_y[2])*par('usr')[3] + (legend_y[2])*par('usr')[4]
    if( diff(legend_y) > diff(legend_x) ){
      align = c("lt","rb")[2]
      gradient = c("x","y")[2]
    }else{
      align = c("lt","rb")[1]
      gradient = c("x","y")[1]
    }  
    
    plotrix::color.legend(xl=xl, yb=yb, xr=xr, yt=yt, legend=round(seq(min(breakpoints),max(breakpoints),length=4),2),
                          rect.col=col(length(breakpoints)-1), cex=0.6, align=align, gradient=gradient)
    
    
    
  }
}
dev.off()
