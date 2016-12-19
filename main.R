library(MASS)
library(mvtnorm)
library(R.matlab)
# data <- readMat('bibtex_new.mat')
data <- readMat('yeast.mat')
source("R/fma.R")
source("R/onlinefma.R")
source("R/fma.em.alg.R")
source("R/mapClass.R")
source("R/misc.R")
source("R/pi.greco.grad.R")
source("R/pi.greco.hess.R")
y <- data[[2]]
y <- as.matrix(y)
fit = onlinefma(y,k=2,r=1,eps=0.0001,scaling=FALSE)

cov <- data[[1]]
cov <- as.matrix(cov)
fit2 = onlinefma(y,k=2,r=1,x.z=cov,x.w=cov,eps=0.0001,scaling=FALSE)