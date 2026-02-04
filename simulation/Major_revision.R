
# One UOT matching example ------------------------------------------------

source("R/thresh_ot_func.R")
library(reticulate)
use_python("~/Library/r-miniconda-arm64/bin/python", required = T) # choose your own python path
source_python("./python/sinkhorn_unbalanced_tv.py")
# Normal case -------------------------------------------------------------

# create 2D normal array
d=2
n_sim = 100
n0 = n_sim
n1 = n_sim
n = n0+n1
idx0 = sort(sample(1:n,n0))
idx = rep(1,n)
idx[idx0] = 0
idx1 = which(idx == 1)


# ----------- generate sample  ----------- #

# set input normal marginals
a = c(-0.4,0.1) # mean 0
b = c(0.5,-0.3) # mean 1

A = matrix( c(0.7, -0.5,
              -0.5,  1), nrow=2, ncol=2) # var 0

B = matrix( c(1, -0.34,
              -0.34,  0.7), nrow=2, ncol=2) # var 1

X0 = MASS::mvrnorm(n0, mu=a, Sigma=A)
X1 = MASS::mvrnorm(n1, mu=b, Sigma=B)
X = matrix(NA,nrow = n, ncol = d)
X[idx0,] = X0
X[idx1,] = X1
plot_cov(X,idx0,idx1,main="2d")




# ----------- uniform marginal weight  ----------- #

# uniform measure
a1 = rep(1/n1,n1)
a0 = rep(1/n0,n0)
center0 = apply(X0, 2, mean)
center1 = apply(X1, 2, mean)

# distance matrix
M = matrix(NA,n0,n1)
for(i in 1:n0){
  for(j in 1:n1){
    M[i,j] = sqrt(sum((X0[i,] - X1[j,])^2))
  }
}
M = M/max(M)
M = t(M)

# UOT estimator with uniform measure
ot <- import("ot")
reg_m_kl = 1e-2/5
reg = 1e-2
maxiter = 2000
ubG = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl, div='kl',
                             numItermax=as.integer(maxiter))

otG = ot$sinkhorn(a1, a0, M, reg, numItermax=as.integer(maxiter))

reg_m_kl2 = 1e-3
ubG2 = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl, div='kl',
                             numItermax=as.integer(maxiter))

# ==== 1. Econ example: it's good to have unmatched observations (weight visualization) ======

# for OT
set.seed(3)
reg_m_kl1 = 1e-2
ubG1 = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl1, div='kl',
                             numItermax=as.integer(maxiter))
reg_m_kl2 = 1e-2/5 # rho
ubG2 = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl2, div='kl',
                             numItermax=as.integer(maxiter))

sum(ubG1)
sum(ubG2)
sum(otG)

table(ubG2 > 1e-2)
table(ubG1 > 1e-2)

cex = 1.2
par(mfrow = c(1,3))
top_pair2 = add_top_match(ubG2,X,idx0,idx1,1:3,main=paste("UOT: rho=",reg_m_kl2),cex = cex)
top_pair1 = add_top_match(ubG1,X,idx0,idx1,1:3,main=paste("UOT: rho=",reg_m_kl1),cex = cex)
top_pair_ot = add_top_match(otG,X,idx0,idx1,1:3,main="OT",cex = cex)
par(mfrow = c(1,1))



# ==== 2/1. Simulation: Estimator does not suffer from curse of dim as much compared to knn ======

d = 2
a = c(-0.4,0.1) # mean 0
b = c(0.5,-0.3) # mean 1

A = matrix( c(0.7, -0.5,
              -0.5,  1), nrow=2, ncol=2) # var 0

B = matrix( c(1, -0.34,
              -0.34,  0.7), nrow=2, ncol=2) # var 1

# closed form
true_UOT = true_2Dnorm_UOT_map(a,b,A,B, reg=reg, reg_marg=reg_m_kl)

# simulation
## sample from closed-form
n_sample = 30
X_cf = MASS::mvrnorm(n_sample, mu = true_UOT$mu, Sigma = true_UOT$H)

## sinkhorn
X0 = MASS::mvrnorm(n0, mu=a, Sigma=A)
X1 = MASS::mvrnorm(n1, mu=b, Sigma=B)
X = matrix(NA,nrow = n, ncol = d)
X[idx0,] = X0
X[idx1,] = X1
# uniform measure
a1 = rep(1/n1,n1)
a0 = rep(1/n0,n0)

# distance matrix
M = matrix(NA,n0,n1)
for(i in 1:n0){
  for(j in 1:n1){
    M[i,j] = sqrt(sum((X0[i,] - X1[j,])^2))
  }
}
M = M/max(M)
M = t(M)
# UOT estimator with uniform measure
ot <- import("ot")
reg_m_kl = 1e-2/5
reg = 1e-2
maxiter = 2000
ubG = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl, div='kl',
                             numItermax=as.integer(maxiter))
ubG = ubG/sum(ubG)
# get empirical moments
ubG_empE1 = rep(0, d*2)
ubG_empE2 = matrix(0, 2*d, 2*d)
for(i in 1:n0){
  for(j in 1:n1){
    ubG_empE1[1:d] = X0[i,]*ubG[j,i] + ubG_empE1[1:d]
    ubG_empE1[1:d+d] = X1[j,]*ubG[j,i] + ubG_empE1[1:d+d]
    
    ubG_empE2[1:d,1:d] = X0[i,]%*%t(X0[i,])*ubG[j,i] + ubG_empE2[1:d,1:d]
    ubG_empE2[1:d,1:d+d] = X0[i,]%*%t(X1[j,])*ubG[j,i] + ubG_empE2[1:d,1:d+d]
    ubG_empE2[1:d+d,1:d] = X1[j,]%*%t(X0[i,])*ubG[j,i] + ubG_empE2[1:d+d,1:d]
    ubG_empE2[1:d+d,1:d+d] = X1[j,]%*%t(X1[j,])*ubG[j,i] + ubG_empE2[1:d+d,1:d+d]
  }
}

cbind(cf_E1 = true_UOT$mu, sinkhorn_E1 = ubG_empE1)

cbind(cf_E2 = c(true_UOT$H), sinkhorn_E2 = c(ubG_empE2))



# use MMD -----------------------------------------------------------------

ker = function(x,y){
  exp(-sum((x-y)^2))
}


X_cf = MASS::mvrnorm(n0*n1, mu = true_UOT$mu, Sigma = true_UOT$H)
X_joint = matrix(NA,nrow = n0*n1, ncol = d*2)
ubG_vec = rep(NA,n0*n1)
for(i in 1:n0){
  for(j in 1:n1){
    X_joint[(i-1)*n1+j,] = cbind(t(X0[i,]),t(X1[j,]))
    ubG_vec[(i-1)*n1+j] = ubG[j,i]
  }
}


emp_MMD = function(X,Y,X_weight){
  E_XX = 0
  for(i in 1:nrow(X)){
    for(j in 1:nrow(X)){
      E_XX = E_XX + X_weight[i]*X_weight[j]*ker(X[i,],X[j,])
    }
  }
  E_XY = 0
  for(i in 1:nrow(X)){
    for(j in 1:nrow(Y)){
      E_XY = E_XY + X_weight[i]*ker(X[i,],Y[j,])/nrow(Y)
    }
  }
  E_YY = 0
  for(i in 1:nrow(Y)){
    for(j in 1:nrow(Y)){
      E_YY = E_YY + ker(Y[i,],Y[j,])/nrow(Y)^2
    }
  }
  
  MMD = E_XX - 2*E_XY + E_YY
  return(MMD)
  
}
sum(ubG_vec)
emp_MMD(X_joint, X_cf, ubG_vec)


# ==== 2/2. fix eps, n, vary rho → compare proportion of perfect matches with different radius rho ======

count_match_within_radius = function(G,M,radius,by_quantile = F, measure_lower=1e-3, measure_quantile = 0.1){
  idx_within_radius = which(M<=radius)
  G_within_radius = G[idx_within_radius]
  if(by_quantile){
    G_lower = quantile(G,measure_quantile)
    return(cbind(num_within_radius = length(idx_within_radius),
                 num_match_within_radius = sum(G_within_radius >= G_lower),
                 proportion = sum(G_within_radius >= G_lower)/length(idx_within_radius)))
  }else{
    return(cbind(num_within_radius = length(idx_within_radius),
                 num_match_within_radius = sum(G_within_radius >= measure_lower),
                 proportion = sum(G_within_radius >= measure_lower)/length(idx_within_radius)))
  }
  
}

count_match_within_radius(otG,M, radius = 0.005, measure_lower=1e-2) # reg_m = inf
count_match_within_radius(ubG1,M, radius = 0.3, measure_lower=1e-2) # reg_m = 1e-2
count_match_within_radius(ubG2,M, radius = 0.3, measure_lower=1e-2) # reg_m = 1e-2/5

plot(ubG1,ubG2, xlab = "UOT: reg_m_kl=1e-2", ylab = "UOT: reg_m_kl=1e-2/5");abline(0,1)

# by quantile
# count_match_within_radius(otG,M, radius = 0.005, by_quantile = F, measure_quantile = 0.5) # reg_m = inf
# count_match_within_radius(ubG1,M, radius = 0.002, by_quantile = F, measure_quantile = 0.99) # reg_m = 1e-2
# count_match_within_radius(ubG2,M, radius = 0.001, by_quantile = F, measure_quantile = 0.99) # reg_m = 1e-2/5

# ? the total mass is different
sum(otG)
sum(ubG1)
sum(ubG2)

# previous ----------------------------------------------------------------

cex = 2
par(mfrow=c(1,3))
get_top_pair_range = function(list_of_top_pairs){
  top_pair_rg_all = NULL
  for(i in 1:length(list_of_top_pairs)){
    top_pair = list_of_top_pairs[[i]]
    top_pair_rg = apply(rbind(X[idx1[top_pair[,1]],] ,
                              X[idx0[top_pair[,2]],]),2,function(x){c(min(x),max(x))})
    top_pair_rg_all = cbind(top_pair_rg_all,top_pair_rg)
  }
  return(apply(top_pair_rg_all,1,function(x){c(min(x),max(x))}))
}

all_range = get_top_pair_range(list(top_pair1,top_pair2,top_pair_ot))
marg = 0.5
bool_in_range = apply(X, 1, function(x){(x[1]<=all_range[2,1]+marg & x[1]>=all_range[1,1]-marg &
                                           x[2]<=all_range[2,2]+marg & x[2]>=all_range[1,2]-marg)})

idx_neighbor = (1:n)[bool_in_range]
par(mfrow = c(1,3))
top_pair1 = add_top_match(ubG1,X,idx0,idx1,1:3,idx_neighbor=idx_neighbor,main=paste("UOT: reg_m_kl=",reg_m_kl1),cex = cex)
top_pair2 = add_top_match(ubG2,X,idx0,idx1,1:3,idx_neighbor=idx_neighbor,main=paste("UOT: reg_m_kl=",reg_m_kl2),cex = cex)
top_pair_ot = add_top_match(otG,X,idx0,idx1,1:3,idx_neighbor=idx_neighbor,main="OT",cex = cex)
par(mfrow = c(1,1))

top_pair_rg = apply(rbind(X[idx1[top_pair[,1]],] ,
                              X[idx0[top_pair[,2]],]),2,function(x){c(min(x),max(x))})
# dist_to_top_pair = apply(X, 1, function(x){sum(abs(x-top_pair_center))})
marg = 0.5
bool_in_range = apply(X, 1, function(x){(x[1]<=top_pair_rg[2,1]+marg & x[1]>=top_pair_rg[1,1]-marg &
                                          x[2]<=top_pair_rg[2,2]+marg & x[2]>=top_pair_rg[1,2]-marg)})

idx_neighbor = (1:n)[bool_in_range]
top_pair = add_top_match(ubG,X,idx0,idx1,idx_neighbor=idx_neighbor,m_vec=1:3,main="UOT",cex = cex)

# top_pair = add_top_match(ubG,X,idx0,idx1,idx_neighbor=NULL,m_vec=1:3,main="OT")
# for KNN3
plot_cov_neighbor(X,idx0,idx1,idx_neighbor=idx_neighbor,main="KNN3",cex)
knn3 = get_knn(M,3)

for(i in 1:dim(top_pair)[1]){
  points(X[idx1[top_pair[i,1]],1],X[idx1[top_pair[i,1]],2], pch = as.character(i),cex = cex)
  knn3_pair = knn3$X1_knn[top_pair[i,1],]
  points(X[idx0[knn3_pair[1]],1], X[idx0[knn3_pair[1]],2], pch = as.character(i),cex = cex)
  points(X[idx0[knn3_pair[2]],1], X[idx0[knn3_pair[2]],2], pch = as.character(i),cex = cex)
  points(X[idx0[knn3_pair[3]],1], X[idx0[knn3_pair[3]],2], pch = as.character(i),cex = cex)
  pair = cbind(X[idx1[top_pair[i,1]],],
               X[idx0[knn3_pair[1]],],
               X[idx0[knn3_pair[2]],],
               X[idx0[knn3_pair[3]],])
  lines(pair[1,c(1,2)],pair[2,c(1,2)])
  lines(pair[1,c(1,3)],pair[2,c(1,3)])
  lines(pair[1,c(1,4)],pair[2,c(1,4)])
}
# knn3$X1_knn[top_pair[1,1],]
# knn1$X1_knn[top_pair[1,1]]
# for KNN1
knn1 = get_knn(M,1)
plot_cov_neighbor(X,idx0,idx1,idx_neighbor=idx_neighbor,dim_plot=dim_plot,main="KNN1",cex = cex)
for(i in 1:dim(top_pair)[1]){
  points(X[idx1[top_pair[i,1]],1],X[idx1[top_pair[i,1]],2], pch = as.character(i),cex = cex)
  knn1_pair = knn1$X1_knn[top_pair[i,1]]
  points(X[idx0[knn1_pair],1], X[idx0[knn1_pair],2], pch = as.character(i),cex = cex)
  pair = cbind(X[idx1[top_pair[i,1]],],
               X[idx0[knn1_pair],])
  lines(pair[1,],pair[2,])
}
par(mfrow=c(1,1))