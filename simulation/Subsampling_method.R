source("R/thresh_ot_func.R")
library(reticulate)
use_python("~/Library/r-miniconda-arm64/bin/python", required = T) # choose your own python path
source_python("./python/sinkhorn_unbalanced_tv.py")
# Normal case -------------------------------------------------------------

# create 2D normal array
d=2
n_sim = 10
n0 = 5*n_sim
n1 = 2*n_sim
n = n0+n1
idx0 = sort(sample(1:n,n0))
idx = rep(1,n)
idx[idx0] = 0
idx1 = which(idx == 1)


# ----------- generate sample  ----------- #

# set input normal marginals
a = c(-0.4,0.1) # mean 0
b = c(0.5,-0.3) # mean 1

A = matrix( c(0.12, -0.35,
              -0.35,  1.5), nrow=2, ncol=2) # var 0

B = matrix( c(1, -0.34,
              -0.34,  0.7), nrow=2, ncol=2) # var 1

X0 = MASS::mvrnorm(n0, mu=a, Sigma=A)
X1 = MASS::mvrnorm(n1, mu=b, Sigma=B)
X = matrix(NA,nrow = n, ncol = d)
X[idx0,] = X0
X[idx1,] = X1
plot_cov(X,idx0,idx1,main="2d")
