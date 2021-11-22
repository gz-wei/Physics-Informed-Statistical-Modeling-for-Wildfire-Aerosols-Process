library(tidyverse)
library(expm)
library(Rfast)
library(fields)
library(MCMCpack)
library(ggpubr)
library(magick)

rm(list = ls()) 

### load the required functions and data 
source(here::here("code", "required_functions.R"))
load(here::here("data", "G17.aod.raw.RData"))
load(here::here("data", "G16.aod.raw.RData"))


### plot the raw data 
### GOES-17
p.G17.aod.true <- list()
for (i in 1:30){
  tempt <- data.frame(long=G17.aod.raw[[i]]$long, lat=G17.aod.raw[[i]]$lat, AOD=G17.aod.raw[[i]]$AOD)
  p.G17.aod.true[[i]] <- 
    ggplot(tempt%>%drop_na(), aes(long, lat, fill=AOD)) +
    geom_raster() + 
    scale_fill_viridis_c(option = "D", limits = c(NA,3.5)) +
    theme(
      panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"),
      panel.grid.major = element_line(size = 0.01, linetype = 'solid',colour = "white"), 
      panel.grid.minor = element_line(size = 0.01, linetype = 'solid',colour = "white"),
      plot.title = element_text(size = 20, hjust=0.5),
      axis.text = element_text(size = 16),
      axis.title=element_text(size=20),
      legend.key.size = unit(0.8, "cm"),
      legend.key.width = unit(0.8,"cm"),
      legend.title = element_text(color = "black", size = 18),
      legend.text = element_text(color = "black", size = 18)
    ) +
    geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
    coord_cartesian(xlim=c(-124,-121.6),ylim=c(35.0,37.4)) +
    ggtitle(paste(c("time", as.character(i)), collapse=" "))
}
print(p.G17.aod.true)
### GOES-16
p.G16.aod.true <- list()
for (i in 1:30){
  tempt <- data.frame(long=G16.aod.raw[[i]]$long, lat=G16.aod.raw[[i]]$lat, AOD=G16.aod.raw[[i]]$AOD)
  p.G16.aod.true[[i]] <- 
    ggplot(tempt%>%drop_na(), aes(long, lat, fill=AOD)) +
    geom_raster() + 
    scale_fill_viridis_c(option = "D", limits = c(NA,3.5)) +
    theme(
      panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"),
      panel.grid.major = element_line(size = 0.01, linetype = 'solid',colour = "white"), 
      panel.grid.minor = element_line(size = 0.01, linetype = 'solid',colour = "white"),
      plot.title = element_text(size = 20, hjust=0.5),
      axis.text = element_text(size = 16),
      axis.title=element_text(size=20),
      legend.key.size = unit(0.8, "cm"),
      legend.key.width = unit(0.8,"cm"),
      legend.title = element_text(color = "black", size = 18),
      legend.text = element_text(color = "black", size = 18)
    ) +
    geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
    coord_cartesian(xlim=c(-124,-121.6),ylim=c(35.0,37.4)) +
    ggtitle(paste(c("time", as.character(i)), collapse=" "))
}
print(p.G16.aod.true)


### calculation the matrix F in the state-space model
N <- 10 ## the highest retained frequency
Nr <- 60 ## the number of vetical and horizontal spatial points of raw AOD images
Omega <- Function_Omega(N)
F <- Function_F(Nr, N, Omega) 


### downsampling and preprosessing for raw data
G17.aod.slt <- UDS(G17.aod.raw, 10)
G16.aod.slt <- UDS(G16.aod.raw, 10)
y1 <- list()
y2 <- list()
for (i in 1:20){
  y1[[i]] <- G17.aod.slt[[i]]$AOD
  y2[[i]] <- G17.aod.slt[[i]]$AOD
}
obs.ccl <- obs_ccl2(y1, y2, F)


### load the physical informtion computed from the optical flow method
load(here::here("data", "K.ifm.RData"))
load(here::here("data", "v.a.RData"))
### Calculation the matrix G in the state-space model
print("this step requires 10-20 mins")
G <- expm(G_ad(1/Nr, v.a, K.ifm, Omega))


### fit the proposed model with Gibbs_FFBS_M2 
m0 <- rep(0.1, 2*N^2)
C0 <- diag(0.01, 2*N^2)
N.sample = 100
start_time <- Sys.time()
fit.M2 <- Gibbs_FFBS_M2(obs.ccl, G, m0, C0, N.sample)
end_time <- Sys.time()
print(end_time-start_time)
image(matrix(F%*%fit.M2$m.flt[[10]][1:100], 60, 60))


### plot the filtering results
dat.map <- map_data("county", "california")
p.G1716.aod.flt <- list()
for(i in 1:20){
  tempt <- data.frame(long=G17.aod.slt[[1]]$long, lat=G17.aod.slt[[1]]$lat, AOD = F%*%fit.M2$m.flt[[i+1]][1:N^2])
  p.G1716.aod.flt[[i]] <- 
    ggplot(tempt, aes(long, lat, fill=AOD)) +
    geom_raster() + 
    scale_fill_viridis_c(option = "D", limits = c(NA,3.5)) +
  theme(
    panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"),
    panel.grid.major = element_line(size = 0.01, linetype = 'solid',colour = "white"), 
    panel.grid.minor = element_line(size = 0.01, linetype = 'solid',colour = "white"),
    plot.title = element_text(size = 20, hjust=0.5),
    axis.text = element_text(size = 16),
    axis.title=element_text(size=20),
    legend.key.size = unit(0.8, "cm"),
    legend.key.width = unit(0.8,"cm"),
    legend.title = element_text(color = "black", size = 18),
    legend.text = element_text(color = "black", size = 18)
  ) + 
    geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
    coord_cartesian(xlim=c(-124,-121.6),ylim=c(35.0,37.4)) +
    ggtitle(paste(c("time", as.character(i)), collapse=" "))
}
print(p.G1716.aod.flt)


### plot the bias correction process
p.G1716.aod.bias <- list()
for(i in 1:20){
  tempt <- data.frame(long=G17.aod.slt[[1]]$long, lat=G17.aod.slt[[1]]$lat, bias = F%*%fit.M2$m.flt[[i+1]][401:800])
  p.G1716.aod.bias[[i]] <- 
    ggplot(tempt%>%drop_na(), aes(long, lat, fill=bias)) +
    geom_raster() + 
    scale_fill_viridis_c(option = "B", limits = c(-2.5,2.5)) +
    theme(
      panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"),
      panel.grid.major = element_line(size = 0.01, linetype = 'solid',colour = "white"),
      panel.grid.minor = element_line(size = 0.01, linetype = 'solid',colour = "white"),
      plot.title = element_text(size = 30, hjust=0.5),
      axis.text = element_text(size = 22),
      axis.title=element_text(size=30),
      legend.key.size = unit(1, "cm"),
      legend.key.width = unit(1,"cm"),
      legend.title = element_text(color = "black", size = 25),
      legend.text = element_text(color = "black", size = 25)
      ,legend.position="none"
    ) +
  geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
    coord_cartesian(xlim=c(-124,-121.6),ylim=c(35.0,37.4)) + 
    ggtitle(paste(c("time", as.character(i)), collapse=" "))
}
print(p.G1716.aod.flt)


### plot the prediction results
p.prd <- list()
G.list <- list()
G.list[[1]] <- G
for(i in 2:10){
  G.list[[i]] = G%*%G.list[[i-1]]
}
for (k in 1:10) {
  data <- data.frame(
    long=G17.aod.slt[[1]]$long,
    lat=G17.aod.slt[[1]]$lat,
    AOD = F%*%G.list[[k]]%*%fit.M2$m.flt[[21]][1:N^2]
  )  
  p.prd[[k]] <- 
    ggplot(data, aes(long, lat, fill=AOD)) +
    geom_raster() + 
    scale_fill_viridis_c(option = "D", limits = c(NA,3.5)) +
    theme(
      panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"),
      panel.grid.major = element_line(size = 0.01, linetype = 'solid',colour = "white"), 
      panel.grid.minor = element_line(size = 0.01, linetype = 'solid',colour = "white"),
      plot.title = element_text(size = 20, hjust=0.5),
      axis.text = element_text(size = 16),
      axis.title=element_text(size=20),
      legend.key.size = unit(0.8, "cm"),
      legend.key.width = unit(0.8,"cm"),
      legend.title = element_text(color = "black", size = 18),
      legend.text = element_text(color = "black", size = 18)
    ) +
  geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
    coord_cartesian(xlim=c(-124,-121.6),ylim=c(35.0,37.4)) +
    ggtitle(paste(c("time", as.character(k+20)), collapse=" "))
}  
print(p.prd)

