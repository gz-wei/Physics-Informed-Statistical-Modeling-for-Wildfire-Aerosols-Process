library(ggplot2)
library(tidyr)
library(expm)
library(MCMCpack)
library(Rfast)
library(matrixStats)
library(Matrix)


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
# plot the initial emission signals and the velocity field in Figure 13
# -------------------------------------------------------------------------------------------------
# the initial emission signals (Figure 13 (a))
N_r <- 20
s <- expand.grid(
  x = seq(0, 1-1/N_r, length.out = N_r),
  y = seq(0, 1-1/N_r, length.out = N_r)
)

dat <- matrix(Source(c(0.1, 0.1), c(0.1,0.2), c(0.2, 0.1), c(0.2,0.2), s), N_r, N_r)
dat.plot <- data.frame(s, signal=c(dat))

ggplot(dat.plot, aes(x, y, fill = signal))+
  geom_raster() +
  theme_paper +
  ggtitle("(a)") + 
  scale_fill_viridis_c(option = "D") 

# plot the velocity field (Figure 13 (b))
angle <- rep(pi/4, 400)
speed <- rep(0.015)
velocity <- data.frame(x=s$x,y=s$y,angle=angle,speed=c(speed))
ggplot(velocity, aes(x, y,fill=speed,angle=angle,radius=scales::rescale(speed, c(0.05, 0.06)))) +
  geom_raster() +
  geom_spoke(arrow = arrow(length = unit(0.035, 'inches'))) +
  scale_fill_viridis_c(option = "D", limits = c(NA,NA)) +
  theme_paper + 
  ggtitle("(b)")


# -------------------------------------------------------------------------------------------------
# obtain the simulated data stream 
# -------------------------------------------------------------------------------------------------
N <- 6  # the highest remaining frequency
Omega <- Function_Omega(N)
coef <- Function_coef(dat, N, Omega)
F <- Function_F(N_r, N, Omega)  # the first input is the recovered pixels (usually equals to N_r)
dat_rc <- pmax(data.frame(s, dat = F%*%coef))

# computation of G matrix without numerical integration 
G <- G_nit(Omega)

# visualize the simulated data stream
noise <- 0.1
N_step <- 29
step <- 0.015
G_step <- expm(step*G)
alpha <- coef
data <- F%*%alpha
p.y_sim <- list()
y_sim <- pmax(data + rnorm(nrow(F), 0, noise),0)
for (i in 1:N_step){
  alpha <- G_step%*%alpha
  tempt <- pmax(F%*%alpha + rnorm(nrow(F), 0, noise),0)
  y_sim <- cbind(y_sim, tempt)
}
for (i in 1:N_step){
  tempt <- data.frame(s, source=y_sim[,i])
  p.y_sim[[i]] <- ggplot(tempt, aes(x,y,fill=source)) +
    geom_raster() +
    theme_paper +
    scale_fill_viridis_c(option = "D", limits = c(0,NA))
  print(p.y_sim[[i]])
  Sys.sleep(0.3)
}

# -------------------------------------------------------------------------------------------------
# Estimated and predictive results from the pure data-driven approach in Figure 14
# -------------------------------------------------------------------------------------------------
y = y_sim[,1:20]
N = 6
N.sample = 4000
T <- ncol(y)
m0 <- rnorm(N^2, 0, 1) 
C0 <- diag(0.1, N^2) 

# MCMC for estimation of G and the augmented state \bm{theta}
sigma2_v <- numeric(N.sample)
W <- replicate(N.sample, matrix(0, N^2, N^2))
G <- replicate(N.sample, matrix(0, N^2, N^2))  # create a storage for MCMC samples
alpha = 8
beta = 1
nu <- N^2
Phi <- diag(0.01, N^2)  
sigma2_v[1] <- 0.01
W[,,1] <- riwish(nu, Phi)
G[,,1] <- riwish(nu, Phi) # initialization

start_time <- Sys.time()
for (k in 2:N.sample) {
  V <- diag(sigma2_v[k-1], nrow(y))
  fit.kf <- Kalman_Filter(y, F, G[,,k-1], V, W[,,k-1], m0, C0)  # forward filtering 

  m.flt <- fit.kf$m.flt 
  m.prd <- fit.kf$m.prd 
  C.flt <- fit.kf$C.flt 
  C.prd <- fit.kf$C.prd
  theta <- matrix(0, N^2, T)
  theta[, T] <- t(rmvnorm(1, m.flt[[T+1]], C.flt[[T+1]]))
  for (t in (T-1):1){
    C.pred.iv <- solve(C.prd[[t+1]])
    ht <- m.flt[[t+1]]+C.flt[[t+1]]%*%t(G[,,k-1])%*%C.pred.iv%*%(theta[,t+1]-m.prd[[t+1]])
    Ht <- C.flt[[t+1]]-C.flt[[t+1]]%*%t(G[,,k-1])%*%C.pred.iv%*%G[,,k-1]%*%C.flt[[t+1]]
    theta[, t] <- t(rmvnorm(1, ht, Ht))
  }  # backward sampling

  beta = beta + sum((y-F%*%theta)^2)/2
  alpha = alpha + T*nrow(y)/2
  sigma2_v[k] <- rinvgamma(1, alpha, beta)  # update sigma2v
  
  Theta <- theta[,1:(T-1)]
  Theta_p <- theta[,2:T]
  r = 8
  Theta.svd <- svd(Theta, nu = r, nv = r)
  Sigma = diag(Theta.svd$d[1:r])
  U = Theta.svd$u
  V = Theta.svd$v
  w <- t(rmvnorm(T-1, rep(0, N^2), W[,,k-1]))
  G[,,k] = (Theta_p-w)%*%V%*%solve(Sigma)%*%t(U)  # calculate G 
  
  nu <- nu+T-1
  Phi <- Phi + (Theta_p-G[,,k]%*%Theta)%*%t(Theta_p-G[,,k]%*%Theta)
  W[,,k] <- riwish(nu, Phi)  # update W
}
end_time <- Sys.time()
# print(end_time-start_time)

# obtain the trained matrix G
G.temp <- diag(0,N^2)
for (i in (3/4*N.sample+1):N.sample){
  G.temp <- G[,,i]+G.temp
}
G.fit <- G.temp/(length((3/4*N.sample+1):N.sample))

# plot the estimated results in Figure 14(a) and (b)
data.plot.a <- data.frame(s, signal = F%*%fit.kf$m.flt[[9]])
ggplot(data.plot.a, aes(x,y,fill=signal)) +
  geom_raster() +
  ggtitle("(a)") + 
  scale_fill_viridis_c(option = "D") +
  theme_paper  # plot Figure 14(a)
data.plot.b <- data.frame(s, signal = F%*%fit.kf$m.flt[[15]])
ggplot(data.plot.b, aes(x,y,fill=signal)) +
  geom_raster() +
  ggtitle("(b)") + 
  scale_fill_viridis_c(option = "D") +
  theme_paper  # plot Figure 14(b)

# plot the estimated results in Figure 14(c) and (d)
data.plot.c <- data.frame(s, signal = F%*%G.fit%*%fit.kf$m.flt[[21]])
ggplot(data.plot.c, aes(x,y,fill=signal)) +
  geom_raster() +
  ggtitle("(c)") + 
  scale_fill_viridis_c(option = "D") +
  theme_paper  # plot Figure 14(c)
data.plot.d <- data.frame(s, signal = F%*%G.fit%*%G.fit%*%G.fit%*%fit.kf$m.flt[[21]])
ggplot(data.plot.d, aes(x,y,fill=signal)) +
  geom_raster() +
  ggtitle("(d)") + 
  scale_fill_viridis_c(option = "D") +
  theme_paper  # plot Figure 14(d)


# -------------------------------------------------------------------------------------------------
# The MSEs of the pure data-driven and the proposed physics-informed approaches (Table 1)
# -------------------------------------------------------------------------------------------------
# for the pure data-driven approach (output: the first row of Table 1)
MSE.dr <- numeric(0)
tempt <- F
for (step in 1:10){
  tempt <- tempt%*%G.fit
  pred <- matrix(tempt%*%fit.kf$m.flt[[21]], 20, 20)
  error <- pred - matrix(y_sim[, step+20], 20,20)
  MSE.dr[step] = sum((error)^2)/nrow(y_sim)
} 
print(MSE.dr)
save(MSE.dr, file = here::here("data", "MSE.dr.RData"))

# for the proposed physics-informed approach (output: from the second row to the last row of Table 1)
K <- matrix(rep(0,400),20,20)
col.dif <- rowDiffs(K)
row.dif <- colDiffs(K)
D.K.x <- cbind(col.dif, col.dif[,19])*20
D.K.y <- rbind(row.dif, row.dif[19,])*20
K.ifm <- data.frame(K=c(K), D.K.x=c(D.K.x), D.K.y=c(D.K.y))

delta <- 0.05
for (step_delta in c(0.015, 0.012, 0.018)){
  id <- c(0, 5, -5, 10, -10, 15, -15, 20, -20, 25, -25, 30, -30)
  G <- list()
  for (i in 1:13){
    v.x <- rep(sqrt(2)*cos((45+id[i])/360*2*pi), (1/delta)^2)
    v.y <- rep(sqrt(2)*sin((45+id[i])/360*2*pi), (1/delta)^2)
    v <- data.frame(v.x=v.x, v.y=v.y)
    G[[i]] <- expm(step_delta*G_ad(delta, v, K.ifm, Omega))
  }
  
  N.sample = 10
  MSE.id = list()
  MSE = list()
  for (j in 1:5) {
    for (i in 1:13) {
      G.tempt = G[[i]]
      y = y_sim[,1:20]
      ## set the parameters for KF function
      T <- ncol(y)
      G.A <- rbind( 
        cbind(G.tempt, diag(N^2)),
        cbind(matrix(0,N^2,N^2), diag(N^2))
      )
      F.A <- cbind(F, matrix(0, nrow(F), ncol(F)))
      m0 <- rnorm(2*N^2, 0, 1) 
      C0 <- diag(0.1, 2*N^2) 
      
      ## create a storage for MCMC samples
      sigma2_v <- numeric(N.sample)
      W <- replicate(N.sample, matrix(0, 2*N^2, 2*N^2))
      ## initialization 
      alpha = 10
      beta = 1
      nu <- 4*N^2
      Phi <- diag(0.01, 2*N^2)  
      sigma2_v[1] <- 0.01
      W[,,1] <- riwish(nu, Phi)
      
      ## Forward filtering using FKF package
      for (k in 2:N.sample) {
        ## Forward filtering 
        V <- diag(sigma2_v[k-1], nrow(y))
        fit.kf <- Kalman_Filter(y, F.A, G.A, V, W[,,k-1], m0, C0)
        
        ## Backward sampling
        m.flt <- fit.kf$m.flt 
        m.prd <- fit.kf$m.prd 
        C.flt <- fit.kf$C.flt 
        C.prd <- fit.kf$C.prd
        theta <- matrix(0, 2*N^2, T)
        theta[, T] <- t(rmvnorm(1, m.flt[[T+1]], C.flt[[T+1]]))
        for (t in (T-1):1){
          C.pred.iv <- solve(C.prd[[t+1]])
          ht <- m.flt[[t+1]]+C.flt[[t+1]]%*%t(G.A)%*%C.pred.iv%*%(theta[,t+1]-m.prd[[t+1]])
          Ht <- C.flt[[t+1]]-C.flt[[t+1]]%*%t(G.A)%*%C.pred.iv%*%G.A%*%C.flt[[t+1]]
          theta[, t] <- t(rmvnorm(1, ht, Ht))
        }
        
        ## update sigma2v
        beta = beta + sum((y-F.A%*%theta)^2)/2
        alpha = alpha + T*nrow(y)/2
        sigma2_v[k] <- rinvgamma(1, alpha, beta)
        
        ## update W
        Theta <- theta[,1:(T-1)]
        Theta_p <- theta[,2:T]
        nu <- nu+T-1
        Phi <- Phi + (Theta_p-G.A%*%Theta)%*%t(Theta_p-G.A%*%Theta)
        W[,,k] <- riwish(nu, Phi)
      }
      
      MSE.tempt <- numeric(0)
      tempt <- F
      for (step in 1:10) {
        tempt <- tempt%*%G.tempt
        pred <- matrix(tempt%*%fit.kf$m.flt[[21]][1:36]+F%*%fit.kf$m.flt[[21]][37:72], 20, 20)
        error <- pred - matrix(y_sim[, step+20], 20,20)
        MSE.tempt[step] = round(sum((error)^2)/nrow(y_sim),digits=4)
      }
      MSE.id[[i]] <- MSE.tempt
    }
    MSE[[j]] <-  MSE.id
  }
  
  MSE.PM.id <- paste(c("MSE.PM", as.character(step_delta*1000)), collapse="")
  assign(MSE.PM.id, MSE) 
  MSE.PM <- replicate(13, matrix(0, 5, 10))
  for (i in 1:5){
    for (j in 1:13){
      MSE.PM[i,,j] <- MSE[[i]][[j]]
    }
  }
  MSE.PM.m = matrix(0, 13, 10)
  for (i in 1:13) {
    MSE.PM.m[i,] = round(colMeans(MSE.PM[,,i]),digits = 4)
  }
  print(MSE.PM.m)
  save(MSE.PM.m, file = here::here("data", paste(c(MSE.PM.id, 'RData'), collapse=".")))
}


# -------------------------------------------------------------------------------------------------
# Prediction MSEs at forward times for Figure 15
# -------------------------------------------------------------------------------------------------
load(here::here("data", "MSE.dr.RData"))  # pure data-driven model 
MSE.dr <- data.frame(
  MSE = MSE.dr,
  Model = rep("Data-driven", length(MSE.dr)),
  time = c("21", '22', '23', '24', '25', '26', '27', '28', '29', '30')
) 

load(here::here("data", "MSE.PM12.RData"))
MSE.PM_12_1 <- MSE.PM.m[1,]  # PI model (0.012, 45)
MSE.PM_12_1 <- data.frame(
  MSE = MSE.PM_12_1,
  Model = rep("(0.012,45)", length(MSE.PM_12_1)),
  time = c("21", '22', '23', '24', '25', '26', '27', '28', '29', '30')
) 

load(here::here("data", "MSE.PM18.RData"))
MSE.PM_18_1 <- MSE.PM.m[1,]  # PI model (0.018, 45)
MSE.PM_18_1 <- data.frame(
  MSE = MSE.PM_18_1,
  Model = rep("(0.018,45)", length(MSE.PM_18_1)),
  time = c("21", '22', '23', '24', '25', '26', '27', '28', '29', '30')
) 

load(here::here("data", "MSE.PM15.RData"))
MSE.PM_15_1 <- MSE.PM.m[1,]  # PI model (0.015, 45)
MSE.PM_15_1 <- data.frame(
  MSE = MSE.PM_15_1,
  Model = rep("(0.015,45)", length(MSE.PM_15_1)),
  time = c("21", '22', '23', '24', '25', '26', '27', '28', '29', '30')
)

MSE.PM_15_6 <- MSE.PM.m[6,]  # PI model (0.015, 60)
MSE.PM_15_6 <- data.frame(
  MSE = MSE.PM_15_6,
  Model = rep("(0.015,60)", length(MSE.PM_15_6)),
  time = c("21", '22', '23', '24', '25', '26', '27', '28', '29', '30')
)

MSE.PM_15_8 <- MSE.PM.m[8,]  # PI model (0.015, 65)
MSE.PM_15_8 <- data.frame(
  MSE = MSE.PM_15_8,
  Model = rep("(0.015,65)", length(MSE.PM_15_8)),
  time = c("21", '22', '23', '24', '25', '26', '27', '28', '29', '30')
)

data.plt <- rbind(MSE.PM_12_1, MSE.PM_15_1, MSE.PM_18_1, MSE.PM_15_6, MSE.PM_15_8, MSE.dr)

# plot Figure 15 (a)
dataset.sub1 <- data.plt %>% 
  subset(time %in% c("21", "22", "23", "24", "25"))
ggplot(dataset.sub1, aes(x=time, y=MSE, group=Model, shape=Model)) +
  geom_line(aes(linetype=Model)) + 
  geom_point(size = 3, aes(shape = Model)) +
  theme_paper

# plot Figure 15 (b)
dataset.sub2 <- data.plt %>% 
  subset(time %in% c("25", "26", "27", "28", "29", "30"))
ggplot(dataset.sub2, aes(x=time, y=MSE, group=Model, shape=Model)) +
  geom_line(aes(linetype=Model)) + 
  geom_point(size = 3, aes(shape = Model)) +
  theme_paper







