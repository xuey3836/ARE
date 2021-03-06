rm(list = ls())
# setwd("C:/Users/o0/Desktop/ARE/")
setwd("/home/o0/Desktop/ARE/")
source("loglik.R")
source("gencc.R")

##set parameter
r = 500
s=  500
k = 0.05
# eta= log(1-1/(exp(k)))
maf = 0.45
lambda = 1.3
thetav = c(0,1/4,1/2,1)
theta= thetav[3]
Nrep= 1000
result = matrix(0,3,Nrep)
are<- function(matdat,theta,thetai,thetaj){
  eta= optimise(f = function(eta){-clogLik(c(1,1,eta),matdat)},interval = c(-10,10))$minimum
  n = sum(matdat)
  L <- chessianLik(c(1,1,eta),matdat)/n
  thetai=  0
  thetaj = 1
  sigma_theta<- function(L,thetai,thetaj){
    A = L[1,3]%*%solve(L[3,3])%*%L[3,1]-L[1,1]
    B = L[2,3]%*%solve(L[3,3])%*%L[3,2]-L[1,2]
    C = L[2,3]%*%solve(L[3,3])%*%L[3,2]-L[2,2]
    sigmah = A*thetai*thetaj+B*(thetai+thetaj)+C
    return(as.numeric(sigmah))
  }
  
  ##s(theta,eta)
  mu_theta_eta = function(theta,eta,mat,lambda){
    n = sum(mat)
    m = length(eta)
    maf= sqrt(colSums(mat)[3]/n)
    mat1 = gencc(r=5000,s=5000,k= 1/(1+exp(-eta)),maf,lambda,theta)
    l1 = cgradLik(c(1,1,eta),mat1)
    L <- chessianLik(c(1,1,eta),mat)/n
    L1 <- chessianLik(c(1,1,eta),mat1)/10000
    s1 = l1[2] + theta*l1[1]
    s2 =(L[2,3]+theta*L[1,3])%*%solve(L[3,3])
    s = s1 - s2%*%l1[3]
    mu = s/n
    mu_1 = (L1[2,2]+theta*L1[1,2]-s2%*%L1[2,3])
    return(list(mu=as.numeric(mu),mu_1=as.numeric(mu_1)))
  }
  
  
  st00 = sigma_theta(L,theta,theta)
  st0i =sigma_theta(L,thetai,theta)
  st0j= sigma_theta(L,thetaj,theta)
  stii=sigma_theta(L,thetai,thetai)
  stjj=sigma_theta(L,thetaj,thetaj)
  stij= sigma_theta(L,thetai,thetaj)
  ##Pitman ARE
  rho_theta0i= st0i/sqrt(st00*stii)
  rho_theta0j= st0j/sqrt(st00*stjj)
  rho_thetaij= stij/sqrt(stii*stjj)
  Ep_MZ= (rho_theta0i+rho_theta0j)^2/(2*(1+rho_thetaij))
  mu_1i = mu_theta_eta(thetai,eta,matdat,lambda=1)$mu_1
  mu_1j = mu_theta_eta(thetaj,eta,matdat,lambda=1)$mu_1
  
  if (theta == thetai){
    mu_1  = mu_1i
  }else if(theta == thetaj){
    mu_1= mu_1j
  }else{
    mu_1 = mu_theta_eta(theta,eta,matdat,lambda=1)$mu_1
  }
  
  ##Chernoff ARE
  Q_Z = 2*(1- pnorm(mu_1/(2*st00)))
  Q_M =2*(1- pnorm((mu_1i/stii+mu_1j/stjj)/
                     (sqrt(8*(1+rho_thetaij)))))
  Ec_MZ = Q_Z/Q_M
  
  ##Hodges-Ledman ARE
  # mu11= mu_theta_eta(theta,eta,matdat,lambda=lambda)
  # mu1i = mu_theta_eta(thetai,eta,matdat,lambda=lambda)
  # mu1j = mu_theta_eta(thetaj,eta,matdat,lambda=lambda)
  # mu = mu11$mu
  # mu_i= mu1i$mu
  # mu_j = mu1j$mu
  # d_Z = (mu/st00)^2
  # d_M =((mu_i/stii+mu_j/stjj)/sqrt(2*(1+rho_thetaij)))^2
  # Ehl_MZ = d_M/d_Z
  # 
  ##Bahadur ARE
  mu_1M = (mu_1i/stii+  mu_1j/stjj)/sqrt(2*(1+rho_thetaij))
  c_Z = 1- pnorm(mu_1/st00)
  c_M = 1- pnorm(mu_1M)
  Eb_MZ = c_Z/c_M
  return(c(Ep_MZ,Ec_MZ,Eb_MZ))
}
thetai = 0
thetaj = 1
for (t in 1:Nrep){
  matdat = gencc(r,s,k,maf,lambda,theta)
  result[,t]= are(matdat = matdat,theta,thetai,thetaj)
}
e_mean=round(rowMeans(result,na.rm=TRUE),5)
e_sd=round(apply(result,1,sd,na.rm=TRUE),5)
e_median=round(apply(result,1, median,na.rm=TRUE),5)
re = rbind(e_mean,e_sd)
View(re)