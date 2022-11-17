library(ggplot2)
library(tidyverse)
library(expm)
library(MCMCpack)
library(Rfast)
library(SpatialVx)
library(bnstruct)
library(matrixStats)


rm(list = ls()) 

source(here::here("util.R"))
load(here::here("data", "G17.aod.raw.RData"))
load(here::here("data", "G16.aod.raw.RData"))

theme_paper <- theme(
  panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"),
  panel.grid.major = element_line(size = 0.01, linetype = 'solid',colour = "white"), 
  panel.grid.minor = element_line(size = 0.01, linetype = 'solid',colour = "white"),
  plot.title = element_text(size = 18, hjust=0.5),
  axis.text = element_text(size = 14),
  axis.title=element_text(size=18),
  legend.key.size = unit(0.8, "cm"),
  legend.key.width = unit(0.8,"cm"),
  legend.title = element_text(color = "black", size = 18),
  legend.text = element_text(color = "black", size = 18)
)


# ------------------------------------------------------------------------------------------------------------------------------------
# plot the velocity and diffusivity fields in figure 9 (need several mins)
# ------------------------------------------------------------------------------------------------------------------------------------
Nr <- 60
s <- expand.grid(
  x = seq(0, 1-1/Nr, length.out = Nr),
  y = seq(0, 1-1/Nr, length.out = Nr)
)

G16.aod.impute <- replicate(60, matrix(0, Nr, Nr))
for (i in 1:60){
  tempt <- matrix(G16.aod.raw[[i]]$AOD, Nr, Nr)
  G16.aod.impute[,,i] <- knn.impute(
    tempt,
    k = 10,
    cat.var = 1:ncol(tempt),
    to.impute = 1:nrow(tempt),
    using = 1:nrow(tempt)
  )
}

dat.velocity <- list()
N.img <- 20
for (i in 1:N.img){
  initial <- G16.aod.impute[,,i]
  final <- G16.aod.impute[,,i+1]
  of.fit <- OF(final, xhat=initial, W=15, verbose=TRUE)
  speed <- matrix(of.fit$err.mag.lin, Nr, Nr)
  angle <- matrix(of.fit$err.ang.lin, Nr, Nr)
  dat.velocity[[i]] <- data.frame(s, speed=c(speed), angle=c(angle/180*pi))
}

dat.velocity.avg.smooth <- velocity(dat.velocity, 1, 19)

v.a.y <- dat.velocity.avg.smooth$speed*sin(dat.velocity.avg.smooth$angle)
v.a.x <- dat.velocity.avg.smooth$speed*cos(dat.velocity.avg.smooth$angle)
v.a <- data.frame(v.a.x=c(v.a.x), v.a.y=c(v.a.y))

speed <- matrix(dat.velocity.avg.smooth$speed,60,60)[seq(1,60,3), seq(1,60,3)]*60*0.04*111*12
long <- matrix(G16.aod.raw[[1]]$long,60,60)[seq(1,60,3), seq(1,60,3)]
lat <- matrix(G16.aod.raw[[1]]$lat,60,60)[seq(1,60,3), seq(1,60,3)]
angle <- matrix(dat.velocity.avg.smooth$angle,60,60)[seq(1,60,3), seq(1,60,3)]
velocity.plot <- data.frame(long=c(long),lat=c(lat),angle=c(angle),speed=c(speed))

dat.map <- map_data("county", "california")
ggplot(velocity.plot, aes(long, lat,fill=speed,angle=angle,radius=scales::rescale(speed, c(0.05, 0.1)))) +
  geom_raster() +
  geom_spoke(arrow = arrow(length = unit(0.035, 'inches'))) +
  scale_fill_viridis_c(option = "D", limits = c(0,30)) +
  labs(fill = expression(paste('speed: ', km/h))) +
  geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
  coord_cartesian(xlim=c(-124,-121.6),ylim=c(35.0,37.4)) + 
  theme_paper  # plot the velocity field 

speed <- matrix(dat.velocity.avg.smooth$speed,60,60)*60*0.04*111*12
long <- matrix(G16.aod.raw[[1]]$long,60,60)
lat <- matrix(G16.aod.raw[[1]]$lat,60,60)
angle <- matrix(dat.velocity.avg.smooth$angle,60,60)
velocity.a <- data.frame(long=c(long),lat=c(lat),angle=c(angle),speed=c(speed))

v.a.y <- velocity.a$speed*sin(velocity.a$angle)
v.a.x <- velocity.a$speed*cos(velocity.a$angle)
v.a <- data.frame(v.a.x=c(v.a.x), v.a.y=c(v.a.y))
v.a.x = matrix(v.a.x, 60, 60)
v.a.y = matrix(v.a.y, 60, 60)
col.dif.x <- rowDiffs(v.a.x)/(0.04*111)
row.dif.x <- colDiffs(v.a.x)/(0.04*111)
col.dif.y <- rowDiffs(v.a.y)/(0.04*111)
row.dif.y <- colDiffs(v.a.y)/(0.04*111)
p1 <- cbind(col.dif.x, rep(NA,60))
p2 <- rbind(row.dif.y, rep(NA,60))
p3 <- rbind(colDiffs(v.a.x), rep(NA,60))
p4 <- cbind(rowDiffs(v.a.y), rep(NA,60))

K <- 0.28*0.04*111*0.04*111*sqrt((p1-p2)^2+(p3+p4)^2)
K.smooth <- image.smooth(K)$z
dat.K <- data.frame(long=c(long), lat=c(lat), diffusivity=c(K.smooth))
ggplot(dat.K, aes(long, lat,fill=diffusivity)) +
  geom_raster() +
  scale_fill_viridis_c(option = "D", limits = c(0,35)) + 
  labs(fill = expression(paste('diffusivity: ', km^2/h))) +
  geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
  coord_cartesian(xlim=c(-124,-121.6),ylim=c(35.0,37.4)) +
  theme_paper  # plot the diffusivity field


# ------------------------------------------------------------------------------------------------------------------------------------
# plot the observed, filtered AOD, and the bias correction process in Figure 8
# ------------------------------------------------------------------------------------------------------------------------------------
# plot row (a) in Figure 8 (observed AOD for GOES 17)
p.G17.aod <- list()
for(i in 1:20){
  tempt <- data.frame(long=G17.aod.raw[[i]]$long, lat=G17.aod.raw[[i]]$lat, AOD = G17.aod.raw[[i]]$AOD)
  p.G17.aod[[i]] <- 
    ggplot(tempt, aes(long, lat, fill=AOD)) +
    geom_raster() + 
    scale_fill_viridis_c(option = "D", limits = c(NA,3.5)) + 
    theme_paper + 
    geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
    coord_cartesian(xlim=c(-124,-121.6),ylim=c(35.0,37.4)) +
    ggtitle(paste(c("time", as.character(i)), collapse=" "))
  print(p.G17.aod[i])
}  

# plot row (b) in Figure 8 (observed AOD for GOES 16)
p.G16.aod <- list()
for(i in 1:20){
  tempt <- data.frame(long=G16.aod.raw[[i]]$long, lat=G16.aod.raw[[i]]$lat, AOD = G16.aod.raw[[i]]$AOD)
  p.G16.aod[[i]] <- 
    ggplot(tempt, aes(long, lat, fill=AOD)) +
    geom_raster() + 
    scale_fill_viridis_c(option = "D", limits = c(NA,3.5)) + 
    theme_paper + 
    geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
    coord_cartesian(xlim=c(-124,-121.6),ylim=c(35.0,37.4)) +
    ggtitle(paste(c("time", as.character(i)), collapse=" "))
  print(p.G16.aod[i])
}  

# plot the row(c) in Figure 8 (filtered AOD) 
# calculation the matrix F 
N <- 20  # the highest retained frequency
Nr <- 60  # the number of vertical and horizontal spatial points of raw AOD images
Omega <- Function_Omega(N)
F <- Function_F(Nr, N, Omega)  
# down-sampling and pre-processing for raw data
G17.aod.slt <- UDS(G17.aod.raw, 20)
G16.aod.slt <- UDS(G16.aod.raw, 20)
y1 <- list()
y2 <- list()
for (i in 1:20){
  y1[[i]] <- G17.aod.slt[[i]]$AOD
  y2[[i]] <- G16.aod.slt[[i]]$AOD
}
obs.ccl <- OBS_ccl2(y1, y2, F)
# Calculation the matrix G 
# load(here::here("data", "K.ifm.RData"))
# load(here::here("data", "v.a.RData"))  # load the physical information computed from the optical flow method
# print("this step requires several hours")
# start_time <- Sys.time()
# G <- expm(G_ad(1/Nr, v.a, K.ifm, Omega))
# save(G, file=here::here("data", "G.RData"))
# end_time <- Sys.time()
# print(end_time-start_time)
load(here::here("data", "G.RData"))
# fit the proposed model with Gibbs_FFBS2
# m0 <- rep(0.1, 2*N^2)
# C0 <- diag(0.01, 2*N^2)
# N.sample = 10
# print("fit the proposed model using Gibbs sampling with FFBS method (10-20 mins)")
# start_time <- Sys.time()
# fit.M2 <- Gibbs_FFBS2(obs.ccl, G, m0, C0, N.sample)
# save(fit.M2, file=here::here("data", "fit.M2.RData"))
# end_time <- Sys.time()
# print(end_time-start_time)
load(here::here("data", "fit.M2.RData"))

p.G1716.aod.flt <- list()
for(i in 1:20){
  tempt <- data.frame(long=G17.aod.slt[[1]]$long, lat=G17.aod.slt[[1]]$lat, AOD = F%*%fit.M2$m.flt[[i+1]][1:N^2])
  p.G1716.aod.flt[[i]] <- 
    ggplot(tempt, aes(long, lat, fill=AOD)) +
    geom_raster() + 
    scale_fill_viridis_c(option = "D", limits = c(NA,3.5)) + 
    theme_paper + 
    geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
    coord_cartesian(xlim=c(-124,-121.6),ylim=c(35.0,37.4)) +
    ggtitle(paste(c("time", as.character(i)), collapse=" "))
  print(p.G1716.aod.flt[i])
}

# plot the row(c) in Figure 8 (filtered AOD the bias correction process)
p.G1716.aod.bias <- list()
for(i in 1:20){
  tempt <- data.frame(long=G17.aod.slt[[1]]$long, lat=G17.aod.slt[[1]]$lat, bias = F%*%fit.M2$m.flt[[i+1]][(N^2+1):(2*N^2)])
  p.G1716.aod.bias[[i]] <- 
    ggplot(tempt%>%drop_na(), aes(long, lat, fill=bias)) +
    geom_raster() + 
    scale_fill_viridis_c(option = "plasma", limits = c(-2.5, 2.5)) +
    geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
    coord_cartesian(xlim=c(-124,-121.6),ylim=c(35.0,37.4)) + 
    ggtitle(paste(c("time", as.character(i)), collapse=" ")) +
    theme_paper
  print(p.G1716.aod.bias[[i]])
}


# ------------------------------------------------------------------------------------------------------------------------------------
# plot the observed AOD and the predicted difference in Figure 10
# ------------------------------------------------------------------------------------------------------------------------------------
# plot row (a) in Figure 10 (observed AOD for GOES 17 at predictive times)
p.G17.aod <- list()
for(i in 21:30){
  tempt <- data.frame(long=G17.aod.raw[[i]]$long, lat=G17.aod.raw[[i]]$lat, AOD = G17.aod.raw[[i]]$AOD)
  p.G17.aod[[i]] <- 
    ggplot(tempt, aes(long, lat, fill=AOD)) +
    geom_raster() + 
    scale_fill_viridis_c(option = "D", limits = c(NA,3.5)) + 
    theme_paper + 
    geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
    coord_cartesian(xlim=c(-124,-121.6),ylim=c(35.0,37.4)) +
    ggtitle(paste(c("time", as.character(i)), collapse=" "))
  print(p.G17.aod[i])
}  

# plot row (c) in Figure 10 (observed AOD for GOES 16 at predictive times)
p.G16.aod <- list()
for(i in 21:25){
  tempt <- data.frame(long=G16.aod.raw[[i]]$long, lat=G16.aod.raw[[i]]$lat, AOD = G16.aod.raw[[i]]$AOD)
  p.G16.aod[[i]] <- 
    ggplot(tempt, aes(long, lat, fill=AOD)) +
    geom_raster() + 
    scale_fill_viridis_c(option = "D", limits = c(NA,3.5)) + 
    theme_paper + 
    geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
    coord_cartesian(xlim=c(-124,-121.6),ylim=c(35.0,37.4)) +
    ggtitle(paste(c("time", as.character(i)), collapse=" "))
  print(p.G16.aod[i])
}  

# plot row (b) in Figure 10 (predictive difference for GOES 17)
p.G1716.prd <- list()
G.list <- list()
G.list[[1]] <- G
for(i in 2:10){
  G.list[[i]] = G%*%G.list[[i-1]]
}
for (k in 1:5) {
  data <- data.frame(
    long=G17.aod.slt[[1]]$long,
    lat=G17.aod.slt[[1]]$lat,
    AOD = F%*%G.list[[k]]%*%fit.M2$m.flt[[21]][1:N^2] - G17.aod.raw[[20+i]]$AOD
  )  
  p.G1716.prd[[k]] <- 
    ggplot(data, aes(long, lat, fill=AOD)) +
    geom_raster() + 
    scale_fill_viridis_c(option = "plasma", limits = c(-3,3)) +
    geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
    coord_cartesian(xlim=c(-124,-121.6),ylim=c(35.0,37.4)) +
    ggtitle(paste(c("time", as.character(k+20)), collapse=" ")) + 
    theme_paper
  print(p.G1716.prd[[k]])
}  

# plot row (d) in Figure 10 (predictive difference for GOES 16)
for (k in 1:5) {
  data <- data.frame(
    long=G17.aod.slt[[1]]$long,
    lat=G17.aod.slt[[1]]$lat,
    AOD = F%*%G.list[[k]]%*%fit.M2$m.flt[[21]][1:N^2] - G16.aod.raw[[20+i]]$AOD
  )  
  p.G1716.prd[[k]] <- 
    ggplot(data, aes(long, lat, fill=AOD)) +
    geom_raster() + 
    scale_fill_viridis_c(option = "plasma", limits = c(-3,3)) +
    geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
    coord_cartesian(xlim=c(-124,-121.6),ylim=c(35.0,37.4)) +
    ggtitle(paste(c("time", as.character(k+20)), collapse=" ")) + 
    theme_paper
  print(p.G1716.prd[[k]])
}  


# ------------------------------------------------------------------------------------------------------------------------------------
# plot the uniform down-sampling in Figure 11
# ------------------------------------------------------------------------------------------------------------------------------------
ds_plot <- UDS(G17.aod.raw, 20)[[1]]
ggplot(ds_plot, aes(long, lat, fill=AOD)) +
  geom_raster() + 
  scale_fill_viridis_c(option = "D", limits = c(NA,3.5), na.value="white") + 
  theme_paper + 
  geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
  coord_cartesian(xlim=c(-124,-121.6),ylim=c(35.0,37.4))


# ------------------------------------------------------------------------------------------------------------------------------------
# plot the computational time v.s. MSE in Figure 12 (several mins)
# ------------------------------------------------------------------------------------------------------------------------------------
for (n in c(20, 21, 24, 28, 29)){
  G17.aod.slt_tempt <- UDS(G17.aod.raw, n)
  G16.aod.slt_tempt <- UDS(G16.aod.raw, n)
  y1_tempt <- list()
  y2_tempt <- list()
  for (i in 1:20){
    y1_tempt[[i]] <- G17.aod.slt_tempt[[i]]$AOD
    y2_tempt[[i]] <- G17.aod.slt_tempt[[i]]$AOD
  }
  obs.ccl_tempt <- OBS_ccl2(y1_tempt, y2_tempt, F)
  # fit the Gibbs_FFBS2
  m0 <- rep(0.1, 2*N^2)
  C0 <- diag(0.01, 2*N^2)
  N.sample_tempt = 2
  start_time <- Sys.time()
  fit.M2 <- Gibbs_FFBS2(obs.ccl_tempt, G, m0, C0, N.sample_tempt)
  end_time <- Sys.time()
  fit.time <- as.numeric(end_time-start_time, units = "mins")
  print(fit.time)
  
  file.name.fit <- paste(c("fit.M2", as.character(n)), collapse=".")
  file.name.time <- paste(c("fit.M2.time", as.character(n)), collapse=".")
  assign(file.name.fit, fit.M2)  
  assign(file.name.time, fit.time)
}

F.A <- cbind(F, matrix(0, nrow(F), ncol(F)))
ms.id <- which(is.na(G16.aod.raw[[20]]$AOD)*1==1)
G16 <- G16.aod.raw[[20]]$AOD[-ms.id]
G17 <- G17.aod.raw[[20]]$AOD

size_id = c(20, 21, 24, 28, 29)
error <- numeric(5)
for (i in 1:5){
  tempt <- get(paste(c("fit.M2", as.character(size_id[i])), collapse="."))
  error[i] <- (sum((G17-F.A%*%tempt$m.flt[[21]])^2) + sum((G16-(F.A%*%tempt$m.flt[[21]])[-ms.id])^2))/(7200-length(ms.id))
}

dat.plot <- data.frame(
  MSE = error*8,
  Size = c(20^2,21^2,24^2,28^2,29^2),
  Time = c(fit.M2.time.20, fit.M2.time.21, fit.M2.time.24, fit.M2.time.28, fit.M2.time.29)/fit.M2.time.29
)

ggplot(dat.plot, aes(x=Size)) +
  geom_line(aes(y = MSE), size=1, linetype = "dashed")+
  geom_line(aes(y = Time), size=1) +
  geom_point(aes(x = Size, y=Time)) +
  geom_point(aes(x = Size, y=MSE)) +
  scale_y_continuous(
    name = "Relative time",
    sec.axis = sec_axis(~.*1/8, name="MSE")
  ) +
  theme_paper 
