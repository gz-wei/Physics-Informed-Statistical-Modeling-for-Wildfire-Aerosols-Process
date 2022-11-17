library(ggplot2)
library(tidyverse)
library(ncdf4)


rm(list = ls()) 

source(here::here("util.R"))

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


# -------------------------------------------------------------------------------------------------
# plot the motivating example in figure 2 
# -------------------------------------------------------------------------------------------------
dat.map <- map_data("county", "california")
long.range <- range(dat.map$long)
lat.range <- range(dat.map$lat)
resolution <- 0.1

# plot for GOES-16 (figure 2 left)
latG16 <- read.table(here::here('data', 'g16_lats.txt'))
longG16 <- read.table(here::here('data', 'g16_lons.txt'))
file.G16.path = here::here('data', 'G16.nc')
dat.G16.plot = Regrid(long.range, lat.range, latG16, longG16, resolution, file.G16.path)

ggplot(dat.G16.plot%>%drop_na(), aes(long, lat, fill = AOD)) +
  geom_raster() +
  scale_fill_viridis_c(option = "D", limits = c(0,3.5)) +
  geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
  ggtitle("GOES-16 at 2020-10-1-18:03") + 
  theme_paper

# plot for GOES-17 (figure 2 right)
latG17 <- read.table(here::here('data', 'g17_lats.txt'))
longG17 <- read.table(here::here('data', 'g17_lons.txt'))
file.G17.path = here::here('data', 'G17.nc')
dat.G17.plot = Regrid(long.range, lat.range, latG17, longG17, resolution, file.G17.path)

ggplot(dat.G17.plot%>%drop_na(), aes(long, lat, fill = AOD)) +
  geom_raster() +
  scale_fill_viridis_c(option = "D", limits = c(0,3.5)) +
  geom_polygon(data = dat.map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill=NA, color = "black") +
  ggtitle("GOES-17 at 2020-10-1-18:03") + 
  theme_paper


# -------------------------------------------------------------------------------------------------
# illustration of the truncation of Fourier terms in Figure 5 (takes several mins)
# -------------------------------------------------------------------------------------------------
load(here::here("data","G17.aod.raw.RData"))

N_r <- sqrt(nrow(G17.aod.raw[[13]]))
s <- expand.grid(
  x = seq(0, 1-1/N_r, length.out = N_r),
  y = seq(0, 1-1/N_r, length.out = N_r)
)
dat <- matrix(G17.aod.raw[[13]]$AOD,N_r,N_r)
data.plot <- data.frame(s,AOD=G17.aod.raw[[13]]$AOD)

plot.initial <-
  ggplot(data.plot, aes(x,y,fill=AOD)) +
  geom_raster() +
  scale_fill_viridis_c(option = "D", limits = c(NA,3.5)) +
  theme_paper +
  ggtitle('(a)')

plot.IFT <- list()
titles <- c('(b)', '(c)', '(d)')
K <- c(40, 20 , 10)
for (i in c(1, 2, 3)){
  N = K[i]
  Omega <- Function_Omega(N)
  coef <- Function_coef(dat, N, Omega)
  F <- Function_F(N_r, N, Omega) 
  dat_rc <- data.frame(s, AOD = F%*%coef)
  plot.IFT[[i]] <-
    ggplot(dat_rc, aes(x, y, fill = AOD))+
    geom_raster() +
    geom_raster() +
    scale_fill_viridis_c(option = "D", limits = c(NA,3.5)) +
    theme_paper + 
    ggtitle(titles[i])
}
  
plot.initial  # plot Figure 5(a)
plot.IFT[[1]]  # plot Figure 5(b)
plot.IFT[[2]]  # plot Figure 5(c)
plot.IFT[[3]]  # plot Figure 5(d)


# -------------------------------------------------------------------------------------------------
# plot the selected area and interested area in Figure 6
# -------------------------------------------------------------------------------------------------
dat.map.fire <- map_data("county", "california")%>%
  subset(subregion%in%c("napa","sonoma", "mendocino", "lake", "solano"))
glass_fire <- ggplot(dat.G17.plot%>%drop_na(), aes(long, lat, fill = AOD)) +
  geom_raster() +
  scale_fill_viridis_c(option = "D", limits = c(0,5)) +
  geom_polygon(data = dat.map, aes(x = long, y = lat, group = group),
               inherit.aes = FALSE, fill=NA, color = "black") +
  geom_polygon(data=dat.map.fire, aes(x=long, y=lat, group=group), inherit.aes=FALSE, fill=NA, color="white") +
  annotate("rect", xmin = -124, xmax = -121.6, ymin = 35.0, ymax = 37.4, alpha = 0.3, fill = "NA", color = "white") +
  theme_paper 
glass_fire


