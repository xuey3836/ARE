
clogLik <- function(beta,mat){
  rv = mat[2,]
  nv = colSums(mat)
  x1 = beta[1]
  x2 = beta[2]
  x3 = beta[3]
  a = rv[2]*log(x1)+rv[3]*log(x2)+sum(rv)*x3
  b = nv[1]*log(1+exp(x3))+nv[2]*log(1+x1*exp(x3))+nv[3]*log(1+x2*exp(x3))
  l = as.numeric(a)-as.numeric(b)
  return(l)
}
# clogLik(beta = c(exp(0.4963),exp(2*0.4963), -0.4909))
cgradLik= function(beta,mat){
  rv = mat[2,]
  nv = colSums(mat)
  x1 = beta[1]
  x2 = beta[2]
  x3 = beta[3]
  ep3 = exp(x3)
  l1 = as.numeric(rv[2]/x1 - nv[2]*ep3/(1+x1*ep3))
  l2 = as.numeric(rv[3]/x2 - nv[3]*ep3/(1+x2*ep3))
  l3 =as.numeric(sum(rv)- nv[1]*ep3/(1+ep3)-nv[2]*x1*ep3/(1+x1*ep3)-nv[3]*x2*ep3/(1+x2*ep3))
  return(c(l1,l2,l3))
}
chessianLik = function(beta,mat){
  rv = mat[2,]
  nv = colSums(mat)
  x1 = beta[1]
  x2 = beta[2]
  x3 = beta[3]
  ep3 = exp(x3)
  h = matrix(0,3,3)
  h[1,1] = as.numeric(-rv[2]/x1^2 -nv[2]*ep3^2/(1+x1*ep3)^2)##
  h[1,2] = 0
  h[1,3] = as.numeric(- nv[2]*ep3/(1+x1*ep3)^2)
  h[2,2] = as.numeric(-rv[3]/x2^2 - nv[3]*ep3^2/(1+x2*ep3)^2)##something
  h[2,3] = as.numeric(- nv[3]*ep3/(1+x2*ep3)^2)
  h[3,3] = as.numeric(-nv[1]*ep3/(1+ep3)^2-nv[2]*x1*ep3/(1+x1*ep3)^2-nv[3]*x2*ep3/(1+x2*ep3)^2)
  h[2,1] = h[1,2]
  h[3,1] =h[1,3]
  h[3,2] =h[2,3]
  return(h)
}

