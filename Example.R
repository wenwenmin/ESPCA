setwd("~/ESPCA-CODE")

source('fun_SPCA.R')
source('fun_ESPCA.R')
#++++++++++++++++++++++++++++++++++
# The first simulated data
set.seed(10)
u1 = rnorm(100)
set.seed(20)
u2 = rnorm(100)

v1 = c(1,-1,0.7,0.1,-0.5,rep(0,5))
v2 = c(rep(0,5),0.1,-2,-0.5,0.3,0.1)

d1 = 10
d2 = 5

X1 = d1*u1%*%t(v1)
X2 = d2*u2%*%t(v2)

set.seed(30)
X = X1 + X2 + 5*matrix(rnorm(100*10),ncol=10)

edges = list()
edges[[1]] = c(1,2)
edges[[2]] = c(1,3)
edges[[3]] = c(1,5)
edges[[4]] = c(2,3)
edges[[5]] = c(2,4)
edges[[6]] = c(3,4)
edges[[7]] = c(4,5)

edges[[8]] = c(6,7)
edges[[9]] = c(6,10)
edges[[10]] = c(7,8)
edges[[11]] = c(7,10)
edges[[12]] = c(8,9)
edges[[13]] = c(8,10)
edges[[14]] = c(9,10)
#++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++
# PCA
out1 = svd(X, nu = 2, nv = 2)

# SPCA
out2 =  SPCA(X, k = 2, kv = c(5,5))

# ESPCA
out3 =  ESPCA(X, k = 2, edges, k.group=6)

# Result 
PC.dat = cbind(cbind(out1$v, out2$V), out3$V)
colnames(PC.dat) = c("PCA.PC1","PCA.PC2","SPCA.PC1","SPCA.PC2","ESPCA.PC1","ESPCA.PC2")
row.names(PC.dat) = paste("Var",1:10, sep = "")
#++++++++++++++++++++++++++++++++++

print(round(PC.dat,3))