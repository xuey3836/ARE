
gencc = function(r,s,k,maf,lambda,theta){
  lambda1 = lambda
  lambda2 = 1- theta+theta*lambda1
  g = c((1-maf)^2,2*maf*(1-maf),maf^2)
  f = rep(0,3)
  f[3] = k/(lambda2*g[1]+lambda1*g[2]+g[3])
  f[2] = lambda1*f[3]
  f[1]= lambda2*f[3]
  p = g*f/k
  q = g*(1-f)/(1-k)
  rv = rowSums(rmultinom(r,1,prob = p))
  sv = rowSums(rmultinom(s,1,prob = q))
  return(rbind(rv,sv))
}
