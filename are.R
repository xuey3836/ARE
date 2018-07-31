rm(list = ls())
setwd("C:/Users/o0/Desktop/ARE/")
source("loglik.R")
n = 1000
# Y = rbinom(n,1,0.3)
G = sample(c(0,1,2),size =n, prob = c(0.2,0.3,0.5),replace = T)
ex = exp(-0.5+ G*0.5)
p1 = ex/(1+ex)
Y = sapply(1:n, function(i){rbinom(1,1,p1[i])})
l =glm(Y~G,family = binomial(link = "logit"))

# library(maxLik)
# maxLik(logLik = clogLik,grad = cgradLik,hess = chessianLik, start = c(1,1,0.5),Y=Y,G=G)
# maxLik(logLik = clogLik, start = c(1,1,0.3))
# optimise(f = function(eta){-clogLik(c(exp(0.5254),exp(2*0.5254),eta),Y,G)},interval = c(-5,5))

## The asympotic relative effcient of Z_MERT/ Z_theta
eta = rep(1,1)
m = length(eta)

eta = optimise(f = function(eta){-clogLik(c(1,1,eta),Y,G)},interval = c(-5,5))$minimum
n = length(Y)
L <- chessianLik(c(1,1,eta),Y,G)/n
theta = 1/2
thetai= 0
thetaj = 1
sigma_theta<- function(L,thetai,thetaj){
  A = L[2,3:(m+2)]%*%solve(L[3:(m+2),3:(m+2)])%*%L[3:(m+2),2]-L[2,2]
  B = L[1,3:(m+2)]%*%solve(L[3:(m+2),3:(m+2)])%*%L[3:(m+2),1]-L[1,2]
  C = L[1,3:(m+2)]%*%solve(L[3:(m+2),3:(m+2)])%*%L[3:(m+2),1]-L[1,1]
  sigmah = A*thetai*thetaj+B*(thetai+thetaj)+C
  return(sigmah)
}

##s(theta,eta)
mu_theta_eta = function(theta,eta,Y,G){
  n = length(Y)
  m = length(eta)
  l1 = cgradLik(c(1,1,eta),Y,G)
  L <- chessianLik(c(1,1,eta),Y,G)/n
  s1 = l1[1] + theta*l1[2]
  s2 =(L[1,3:(m+2)]+theta*L[2,3:(m+2)])%*%solve(L[3:(m+2),3:(m+2)])
  s = s1 - s2%*%l1[3]
  mu = s/n
  mu_1 = (L[1,1]+theta*L[1,2]-s2%*%L[1,3:(m+2)])/n
  return(list(mu=mu,mu_1=mu_1))
}
##


##Pitman ARE
rho_theta0i= sigma_theta(L,theta,thetai)/(sigma_theta(L,theta,theta)*sigma_theta(L,thetai,thetai))
rho_theta0j= sigma_theta(L,theta,thetaj)/(sigma_theta(L,theta,theta)*sigma_theta(L,thetaj,thetaj))
rho_thetaij= sigma_theta(L,thetai,thetaj)/((sigma_theta(L,thetai,thetai)*sigma_theta(L,thetaj,thetaj)))
Ep_MZ= (rho_theta0i+rho_theta0j)^2/(2*(1+rho_thetaij))

##Chernoff ARE
mu11= mu_theta_eta(theta,eta,Y,G)
mu1i = mu_theta_eta(thetai,eta,Y,G)
mu1j = mu_theta_eta(thetaj,eta,Y,G)
mu = mu11$mu
mu_1 = mu11$mu_1
mu_1i = mu1i$mu_1
mu_ij = mu1j$mu_1
Q_Z = 2(1- pnorm(mu_1/(2*sigma_theta(theta,theta))))
Q_M =2(1- pnorm((mu_1i/(2*sigma_theta(thetai,thetai))+mu_1j/(2*sigma_theta(thetaj,thetaj)))/
                  (sqrt(8*(1+rho_thetaij)))))
Ec_MZ = Q_M/Q_Z

##Hodges-Ledman ARE
d_Z = mu^2/sigma_theta(L,theta,theta)
d_M =( mu_i/sigma_theta(L,thetai,thetai)+  mu_j/sigma_theta(L,thetaj,thetaj))/sqrt(2*(1+rho_thetaij))
E_HL = d_M/d_Z

##Bahadur ARE
mu_1M = ( mu_1i/sigma_theta(L,thetai,thetai)+  mu_1j/sigma_theta(L,thetaj,thetaj))/sqrt(2*(1+rho_thetaij))
c_Z = 1- pnorm(mu_1/sigma_theta(L,theta))
c_M = 1- pnorm(mu_1M)
E_B = c_M/c_Z


