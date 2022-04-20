library(dplyr)

# Load fits ---------------------------------------------------------------

work_dir <- "04_outputs/user_grid_outEBS/no_sex/finescale/c(3,3,3,3,2,2)/no2000/av_season/"

load(paste0(work_dir,"Overlap.RData"))
Overlap_av <- Overlap %>% mutate(Model="averaged")
load(paste0(work_dir,"fit.RData"))
fit_av <-fit

load(paste0("04_outputs/user_grid_outEBS/no_sex/finescale/c(3,3,3,3,2,2)/no2000/sum/Overlap.RData"))
Overlap_sum <- Overlap %>% mutate(Model="sum")
load(paste0("04_outputs/user_grid_outEBS/no_sex/finescale/c(3,3,3,3,2,2)/no2000/sum/fit.RData"))
fit_sum <- fit

Overlap <- bind_rows(Overlap_av,Overlap_sum)

# index1
sd1 <- (summary(fit_av$parameter_estimates[["SD"]])[,2])
sd1 <- as.vector(sd1[which(names(sd1)=="Index_ctl")])

index1 <- as.data.frame(cbind(fit_av$Report$Index_ctl, c(2001:2019),sd1))
index1 <- index1 %>% mutate(Model="averaged")
colnames(index1) <- c("Index", "Years","sd","Model")
index1 <- as_tibble(index1) %>% dplyr::mutate(low=Index-0.67*sd, high=Index+0.67*sd)


# index2
sd2 <- (summary(fit_sum$parameter_estimates[["SD"]])[,2])
sd2 <- as.vector(sd2[which(names(sd2)=="Index_ctl")])

index2 <- as.data.frame(cbind(fit_sum$Report$Index_ctl, c(2001:2019),sd2))
index2 <- index2 %>% mutate(Model="sum")
colnames(index2) <- c("Index", "Years","sd","Model")

index2 <- as_tibble(index2) %>% dplyr::mutate(low=Index-0.67*sd, high=Index+0.67*sd)

Index <- bind_rows(index1, index2)

q <- ggplot() +
  # accrross time
  geom_line(data=Index,aes(x=Years,y=Index,color=Model,linetype=Model),size=2)+
  geom_ribbon(data=Index,aes(x=Years,ymin=low,ymax=high,fill=Model),alpha=0.2)
q

ggsave(paste0(work_dir,'Comp_sum_av.png'), q)


q <- ggplot() +
  # accrross time
  geom_line(data=Overlap,aes(x=Years,y=Index,color=Model,linetype=Model),size=2)

q

ggsave(paste0(work_dir,'Overalpsum_av.png'), q)

