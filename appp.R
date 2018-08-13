library(xlsx)
dat = read.xlsx("application.xls",1)
res = matrix(0,20,3)
for (i in c(2:9,11:16,18:20)){
  # i=3
  datmat = rbind(as.numeric(dat[i,5:7]),dat[i,2:4])
  colnames(datmat)<- c("g0","g1","g2")
  rownames(datmat)= c("sv","rv")
  res[i,] = round(are(matdat = as.matrix(datmat),theta= 1/2, thetai= 0,thetaj=1),4)
  }