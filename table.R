setwd("C:/Users/o0/Desktop/ARE/ref/")
library(xlsx)
a= read.xlsx(file = "result-are0810.xlsx",2)
b= a[12,c(4:6,8:10,12:14)]
c= NULL
for (i in 1:9){
   c0 = paste0(round(as.numeric(as.character(b[1,i])),3))
   if (i%%3==1){
     c = paste0(c,"&&",c0)
   }else{
     c = paste0(c,"&",c0)
   }
}
b1= a[13,c(4:6,8:10,12:14)]
c1= NULL
for (i in 1:9){
  c0 = paste0(round(as.numeric(as.character(b1[1,i])),3))
  if (i%%3==1){
    c1 = paste0(c1,"&&(",c0,")")
  }else{
    c1 = paste0(c1,"&(",c0,")")
  }
}

