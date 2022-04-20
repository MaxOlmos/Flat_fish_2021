# Comparison --------------------------------------------------------------
# -------------------------------------------------------------------------
# male vs female
library(dplyr)

# Load fits ---------------------------------------------------------------

work_dir <- "04_outputs/user_grid_outEBS/sex/male/no2000/"
load(paste0(work_dir,"fit.RData"))
fit_male <-fit
fit$Report

work_dir <- "04_outputs/user_grid_outEBS/sex/female/no2000/"
load(paste0(work_dir,"fit.RData"))
fit_female <-fit

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
  geom_ribbon(data=Index,aes(x=Years,ymin=low,ymax=high,fill=sex),alpha=0.2)+ theme_bw()
q

ggsave(paste0(work_dir,'Comp_male_female_av_season.png'), q)
