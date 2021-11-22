### all required functions 

###Name: Uniformly downsampling for raw dataset
###Author: Guanzhou Wei
###Last revision (dd/mm/yyyy): 08/06/2021
###Required packages: "tidyverse"
UDS <- function(datset, n.slt){
  test <- datset[[1]]
  long.range <- range(test$long)
  lat.range <- range(test$lat)
  lat.slt = seq(lat.range[1],lat.range[2], 0.04)[60/n.slt*(1:n.slt)]
  long.slt = seq(long.range[1],long.range[2], 0.04)[60/n.slt*(1:n.slt)]
  for (i in 1:length(datset)){
    datset[[i]]$id <- 1:nrow(datset[[i]])
    s.slt <- datset[[i]] %>%
      subset(lat%in%round(lat.slt,digits=2) & long%in%round(long.slt,digits=2))
    datset[[i]]$AOD <- NaN
    datset[[i]]$AOD[s.slt$id] <- s.slt$AOD
  }
  return(datset)
}


###Name: Generate the wavenumbe set
Function_Omega <- function(N){
  Omega1 <- data.frame(
    k1 = c(0, 0, N/2, N/2),
    k2 = c(0, N/2, 0, N/2)
  )
  Omega2_p1 <- data.frame(
    k1 = seq(1, N/2-1, 1),
    k2 = N/2
  )
  Omega2_p2 <- expand.grid(
    k1 = seq(0, N/2, 1),
    k2 = seq(1, N/2-1, 1)
  ) 
  Omega2_p3 <- expand.grid(
    k1 = seq(1, N/2-1, 1),
    k2 = seq(-N/2+1, 0, 1)
  ) 
  Omega2 <- rbind(Omega2_p1, Omega2_p2, Omega2_p3)
  return(list(Omega1=Omega1, Omega2=Omega2))
}


###Name: Calculation the matrix F in the state-space model
Function_F <- function(N_rcr, N, Omega){
  Omega1 <- Omega$Omega1
  Omega2 <- Omega$Omega2
  sites_rc <- expand.grid(
    seq(0, 1-1/N_rcr, length.out = N_rcr),
    seq(0, 1-1/N_rcr, length.out = N_rcr)
  )
  F1 <- matrix(0, nrow(sites_rc), nrow(Omega1))
  for (i in 1:nrow(sites_rc)){
    for(j in 1:nrow(Omega1)){
      F1[i,j] <- cos(2*pi*sites_rc[i,1]*Omega1[j,1]+2*pi*sites_rc[i,2]*Omega1[j,2])
    }
  }
  F2c <- matrix(0, nrow(sites_rc), nrow(Omega2))
  for (i in 1:nrow(sites_rc)){
    for(j in 1:nrow(Omega2)){
      F2c[i,j] <- 2*cos(2*pi*sites_rc[i,1]*Omega2[j,1]+2*pi*sites_rc[i,2]*Omega2[j,2])
    }
  }
  F2s <- matrix(0, nrow(sites_rc), nrow(Omega2))
  for (i in 1:nrow(sites_rc)){
    for(j in 1:nrow(Omega2)){
      F2s[i,j] <- 2*sin(2*pi*sites_rc[i,1]*Omega2[j,1]+2*pi*sites_rc[i,2]*Omega2[j,2])
    }
  }
  rv <- cbind(F1, F2c, F2s)
  return(rv)
}


###Name: Trainsition matrix (G) calculation with physical parameters
###Author: Guanzhou Wei
###Last revision (dd/mm/yyyy): 06/06/2021
###Required functions: "Phi1.int.R","Phi2.int.R","Phi3.int.R","Phi5.int.R","Phi6.int.R","Phi7.int.R","Phi8.int.R"
Phi1.int <- function(k1, k2, m1, m2, delta, v) {
  f <- function(s) {2*pi*(s[3]*k1+s[4]*k2)*sin(2*pi*s[1]*k1+2*pi*s[2]*k2)*cos(2*pi*s[1]*m1+2*pi*s[2]*m2)}
  N.row <- 1/delta
  x <- seq(0, 1-1/N.row, length.out = N.row)
  y <- seq(0, 1-1/N.row, length.out = N.row)
  z <- as.matrix(cbind(expand.grid(x, y), v))
  rv <- sum(apply(z, 1, f))*delta^2
  return(rv)
}
Phi2.int <- function(k1, k2, m1, m2, delta, v) {
  f <- function(s) {-2*pi*(s[3]*k1+s[4]*k2)*cos(2*pi*s[1]*k1+2*pi*s[2]*k2)*cos(2*pi*s[1]*m1+2*pi*s[2]*m2)}
  N.row <- 1/delta
  x <- seq(0, 1-1/N.row, length.out = N.row)
  y <- seq(0, 1-1/N.row, length.out = N.row) 
  z <- as.matrix(cbind(expand.grid(x, y), v))
  rv <- sum(apply(z, 1, f))*delta^2
  return(rv)
}
Phi3.int <- function(k1, k2, m1, m2, delta, v) {
  f <- function(s) {2*pi*(s[3]*k1+s[4]*k2)*sin(2*pi*s[1]*k1+2*pi*s[2]*k2)*sin(2*pi*s[1]*m1+2*pi*s[2]*m2)}
  N.row <- 1/delta
  x <- seq(0, 1-1/N.row, length.out = N.row)
  y <- seq(0, 1-1/N.row, length.out = N.row)
  z <- as.matrix(cbind(expand.grid(x, y), v))
  rv <- sum(apply(z, 1, f))*delta^2
  return(rv)
}
Phi5.int <- function(k1, k2, m1, m2, delta, K.ifm) {
  f <- function(s) {(-4*pi^2*(k1^2+k2^2)*s[3]*cos(2*pi*s[1]*k1+2*pi*s[2]*k2)-2*pi*(s[4]*k1+s[5]*k2)*sin(2*pi*s[1]*k1+2*pi*s[2]*k2))*cos(2*pi*s[1]*m1+2*pi*s[2]*m2)}
  N.row <- 1/delta
  x <- seq(0, 1-1/N.row, length.out = N.row)
  y <- seq(0, 1-1/N.row, length.out = N.row)
  z <- as.matrix(cbind(expand.grid(x, y), K.ifm))
  rv <- sum(apply(z, 1, f))*delta^2
  return(rv)
}
Phi6.int <- function(k1, k2, m1, m2, delta, K.ifm) {
  f <- function(s) {(-4*pi^2*(k1^2+k2^2)*s[3]*sin(2*pi*s[1]*k1+2*pi*s[2]*k2)-2*pi*(s[4]*k1+s[5]*k2)*cos(2*pi*s[1]*k1+2*pi*s[2]*k2))*cos(2*pi*s[1]*m1+2*pi*s[2]*m2)}
  N.row <- 1/delta
  x <- seq(0, 1-1/N.row, length.out = N.row)
  y <- seq(0, 1-1/N.row, length.out = N.row)
  z <- as.matrix(cbind(expand.grid(x, y), K.ifm))
  rv <- sum(apply(z, 1, f))*delta^2
  return(rv)
}
Phi7.int <- function(k1, k2, m1, m2, delta, K.ifm) {
  f <- function(s) {(-4*pi^2*(k1^2+k2^2)*s[3]*cos(2*pi*s[1]*k1+2*pi*s[2]*k2)-2*pi*(s[4]*k1+s[5]*k2)*sin(2*pi*s[1]*k1+2*pi*s[2]*k2))*sin(2*pi*s[1]*m1+2*pi*s[2]*m2)}
  N.row <- 1/delta
  x <- seq(0, 1-1/N.row, length.out = N.row)
  y <- seq(0, 1-1/N.row, length.out = N.row)
  z <- as.matrix(cbind(expand.grid(x, y), K.ifm))
  rv <- sum(apply(z, 1, f))*delta^2
  return(rv)
}
Phi8.int <- function(k1, k2, m1, m2, delta, K.ifm) {
  f <- function(s) {(-4*pi^2*(k1^2+k2^2)*s[3]*sin(2*pi*s[1]*k1+2*pi*s[2]*k2)-2*pi*(s[4]*k1+s[5]*k2)*cos(2*pi*s[1]*k1+2*pi*s[2]*k2))*sin(2*pi*s[1]*m1+2*pi*s[2]*m2)}
  N.row <- 1/delta
  x <- seq(0, 1-1/N.row, length.out = N.row)
  y <- seq(0, 1-1/N.row, length.out = N.row)
  z <- as.matrix(cbind(expand.grid(x, y), K.ifm))
  rv <- sum(apply(z, 1, f))*delta^2
  return(rv)
}

G_ad <- function(delta, v, K.ifm, Omega){
  Omega1 <- Omega$Omega1
  Omega2 <- Omega$Omega2
  N_Gr <- 2*nrow(Omega2) + nrow(Omega1)
  G <- matrix(0, N_Gr, N_Gr)
  k1 <- nrow(Omega1)
  k2 <- nrow(Omega2)
  ##for real part in Omega1 
  ##for special case, m = 1
  for (k in 1:k1){
    G[1, k] <- Phi1.int(Omega1[k,1],Omega1[k,2],Omega1[1,1],Omega1[1,2],delta,v) +
      Phi5.int(Omega1[k,1],Omega1[k,2],Omega1[1,1],Omega1[1,2],delta,K.ifm)
  }
  for (k in (k1+1):(k1+k2)){
    G[1, k] <- 2*Phi1.int(Omega2[k-k1,1],Omega2[k-k1,2],Omega1[1,1],Omega1[1,2],delta,v) + 
      2*Phi5.int(Omega2[k-k1,1],Omega2[k-k1,2],Omega1[1,1],Omega1[1,2],delta,K.ifm)
  }
  for (k in (k1+k2+1):(k1+2*k2)){
    G[1, k] <- 2*Phi2.int(Omega2[k-(k1+k2),1],Omega2[k-(k1+k2),2],Omega1[1,1], Omega1[1,2],delta,v) + 
      2*Phi6.int(Omega2[k-(k1+k2),1],Omega2[k-(k1+k2),2],Omega1[1,1], Omega1[1,2],delta,K.ifm)
  }
  ##for m!=1 in Omega1 
  for (m in 2:k1){
    for (k in 1:k1){
      G[m, k] <- 2*Phi1.int(Omega1[k,1],Omega1[k,2],Omega1[m,1],Omega1[m,2],delta,v) + 
        2*Phi5.int(Omega1[k,1],Omega1[k,2],Omega1[m,1],Omega1[m,2],delta,K.ifm)
    }
    for (k in (k1+1):(k1+k2)){
      G[m, k] <- 4*Phi1.int(Omega2[k-k1,1],Omega2[k-k1,2],Omega1[m,1],Omega1[m,2],delta,v) + 
        4*Phi5.int(Omega2[k-k1,1],Omega2[k-k1,2],Omega1[m,1],Omega1[m,2],delta,K.ifm)
    }
    for (k in (k1+k2+1):(k1+2*k2)){
      G[m, k] <- 4*Phi2.int(Omega2[k-(k1+k2),1],Omega2[k-(k1+k2),2],Omega1[m,1],Omega1[m,2],delta,v) + 
        4*Phi6.int(Omega2[k-(k1+k2),1],Omega2[k-(k1+k2),2],Omega1[m,1],Omega1[m,2],delta,K.ifm)
    }
  }
  ##for real part in Omega2 
  for (m in 1:k2){
    for (k in 1:k1){
      G[m+k1, k] <- Phi1.int(Omega1[k,1],Omega1[k,2],Omega2[m,1],Omega2[m,2],delta,v) + 
        Phi5.int(Omega1[k,1],Omega1[k,2],Omega2[m,1],Omega2[m,2],delta,K.ifm)
    }
    for (k in (k1+1):(k1+k2)){
      G[m+k1, k] <- 2*Phi1.int(Omega2[k-k1,1],Omega2[k-k1,2],Omega2[m,1],Omega2[m,2],delta,v) + 
        2*Phi5.int(Omega2[k-k1,1],Omega2[k-k1,2],Omega2[m,1],Omega2[m,2],delta,K.ifm)
    }
    for (k in (k1+k2+1):(k1+2*k2)){
      G[m+k1, k] <- 2*Phi2.int(Omega2[k-(k1+k2),1],Omega2[k-(k1+k2),2],Omega2[m,1],Omega2[m,2],delta,v) +
        2*Phi6.int(Omega2[k-(k1+k2),1],Omega2[k-(k1+k2),2],Omega2[m,1],Omega2[m,2],delta,K.ifm)
    }
  }
  ##for imaginary part in Omega2
  for (m in 1:k2){
    for (k in 1:k1){
      G[m+k1+k2, k] <- Phi3.int(Omega1[k,1],Omega1[k,2],Omega2[m,1],Omega2[m,2],delta,v) + 
        Phi7.int(Omega1[k,1],Omega1[k,2],Omega2[m,1],Omega2[m,2],delta,K.ifm)
    }
    for (k in (k1+1):(k1+k2)){
      G[m+k1+k2, k] <- 2*Phi3.int(Omega2[k-k1,1],Omega2[k-k1,2],Omega2[m,1],Omega2[m,2],delta,v) + 
        2*Phi7.int(Omega2[k-k1,1],Omega2[k-k1,2],Omega2[m,1],Omega2[m,2],delta,K.ifm)
    }
    for (k in (k1+k2+1):(k1+2*k2)){
      G[m+k1+k2, k] <- -2*Phi1.int(Omega2[m,1],Omega2[m,2],Omega2[k-(k1+k2),1],Omega2[k-(k1+k2),2],delta,v) +
        2*Phi8.int(Omega2[m,1],Omega2[m,2],Omega2[k-(k1+k2),1],Omega2[k-(k1+k2),2],delta,K.ifm)
    }
  }
  return(G)
}



###Name: Boolean matrix for missing points
M_id <- function(y){
  ms.id <- which(is.na(y)*1==0)
  M <- matrix(0, length(ms.id), length(y))
  for (i in 1:length(ms.id)){
    M[i, ms.id[i]] <- 1
  }
  return(list(M=M, ms.id=ms.id))
}


###Name: preprosessing the raw data 
###Required functions: "M_id.R"
obs_ccl2 <- function(y1, y2, F){
  T <- length(y1)
  F1t <- list()
  F2t <- list()
  y1.c <- list() 
  y2.c <- list()
  id1 <- list()
  id2 <- list()
  for (i in 1:T){
    tempt1 <- M_id(y1[[i]])
    tempt2 <- M_id(y2[[i]])
    F1t[[i]] <- tempt1$M%*%F
    F2t[[i]] <- tempt2$M%*%F
    y1.c[[i]] <- y1[[i]][tempt1$ms.id]
    y2.c[[i]] <- y2[[i]][tempt2$ms.id]
    id1[[i]] <- tempt1$ms.id
    id2[[i]] <- tempt2$ms.id
  }
  return(list(F1t=F1t,F2t=F2t,y1.c=y1.c,y2.c=y2.c,id1=id1,id2=id2))
}


# Name: Gibbs sampling with FFBS algorithm for parameters estimation of the state-space model
Gibbs_FFBS_M2 <- function(obs.ccl, G, m0, C0, N.sample){
  ### input data 
  y1 <- obs.ccl$y1.c
  y2 <- obs.ccl$y2.c
  T <- length(y1)
  F1t <- obs.ccl$F1t
  F2t <- obs.ccl$F2t
  id1 <- obs.ccl$id1
  id2 <- obs.ccl$id2
  G <- rbind( 
    cbind(G, diag(nrow(G))),
    cbind(matrix(0,nrow(G),nrow(G)), diag(nrow(G)))
  )
  ### initialization and storage for Gibbs sampler 
  sigma12_v <- numeric(N.sample)
  sigma22_v <- numeric(N.sample)
  W <- list()
  alpha1 =  alpha2 = 1
  beta1 = beta2 = 0.1
  nu <- nrow(C0)
  Phi <- diag(0.01, nrow(C0))  
  sigma12_v[1] <- rinvgamma(1, alpha1, beta1)
  sigma22_v[1] <- rinvgamma(1, alpha2, beta2)
  W[[1]] <- rwish(nu, Phi)
  #### loop for Gibbs sampling
  for (k in 2:N.sample){
    ### Forward Filtering
    m.flt = m.prd = C.flt = C.prd = vector("list")
    m.flt[[1]] <- m0
    C.flt[[1]] <- C0
    for (t in 1:T){
      m.prd[[t]] <- G%*%m.flt[[t]]
      C.prd[[t]] <- G%*%C.flt[[t]]%*%t(G) + W[[k-1]]
      sigma2t_v <- c(rep(sigma12_v[k-1], length(y1[[t]])), rep(sigma22_v[k-1], length(y2[[t]])))
      Vt <- diag(sigma2t_v)
      Ft <- rbind(F1t[[t]], F2t[[t]])
      Ftt <- cbind(Ft, matrix(0,nrow(Ft),ncol(Ft)))
      Q <-  Ftt%*%C.prd[[t]]%*%t(Ftt) + Vt 
      Q.inv <- solve(Q)
      yt <- c(y1[[t]], y2[[t]])
      m.flt[[t+1]] <- m.prd[[t]]+C.prd[[t]]%*%t(Ftt)%*%Q.inv%*%(yt-Ftt%*%m.prd[[t]])
      C.flt[[t+1]] <- C.prd[[t]]-C.prd[[t]]%*%t(Ftt)%*%Q.inv%*%Ftt%*%C.prd[[t]]
    }
    m.prd[[T+1]] <- G%*%m.flt[[t]]
    C.prd[[T+1]] <- G%*%C.flt[[t]]%*%t(G) + W[[k-1]]
    ### image(matrix(F%*%(G%*%m.flt[[8]])[1:100], 20, 20))
    ### Backward Sampling
    theta <- matrix(0, nrow(G), T)
    theta[, T] <- t(rmvnorm(1, m.flt[[T+1]], C.flt[[T+1]]))
    for (t in (T-1):1){
      ht <- m.flt[[t+1]]+C.flt[[t+1]]%*%t(G)%*%solve(C.prd[[t+1]])%*%(theta[,t+1]-m.prd[[t+1]])
      Ht <- C.flt[[t+1]]-C.flt[[t+1]]%*%t(G)%*%solve(C.prd[[t+1]])%*%G%*%C.flt[[t+1]]
      theta[, t] <- t(rmvnorm(1, ht, Ht))
    }
    #image(matrix(F%*%(G%*%theta[,20])[1:100], 20, 20))
    ### resudual calculation
    y1.rsd.clss <- rep(0,T)
    y1.rsd.cln <- rep(0,T)
    for (t in 1:T){
      tempt <- y1[[t]]-F1t[[t]]%*%theta[1:(1/2*length(m0)),t]
      y1.rsd.cln[t] <- length(tempt)
      y1.rsd.clss[t] <- sum(tempt^2)
    }
    y2.rsd.clss <- rep(0,T)
    y2.rsd.cln <- rep(0,T)
    for (t in 1:T){
      tempt <- y2[[t]]-F2t[[t]]%*%theta[1:(1/2*length(m0)),t]
      y2.rsd.cln[t] <- length(tempt)
      y2.rsd.clss[t] <- sum(tempt^2)
    }
    state.rsd <- matrix(0, length(m0), T-1)
    for (t in 1:(T-1)){
      state.rsd[,t] <- theta[,t+1] - G%*%theta[,t]
    }
    #### update parameters
    alpha1 <- alpha1 + sum(y1.rsd.cln)/2
    beta1 <- beta1 + sum(y1.rsd.clss)/2
    alpha2 <- alpha2 + sum(y2.rsd.cln)/2
    beta2 <- beta2 + sum(y2.rsd.clss)/2
    nu <- nu + T-1
    Phi <- Phi + state.rsd%*%t(state.rsd)
    sigma12_v[k] <- rinvgamma(1, alpha1, beta1)
    sigma22_v[k] <- rinvgamma(1, alpha2, beta2)
    W[[k]] <- riwish(nu, Phi)
  }
  return(list(W=W, sigma12_v=sigma12_v, sigma22_v=sigma22_v, m.flt=m.flt))
}
