## set input normal marginals ##
n0 = 5*n_sim
n1 = 2*n_sim
n = n0+n1
idx0 = sort(sample(1:n,n0))
idx = rep(1,n)
idx[idx0] = 0
idx1 = which(idx == 1)

X0 = MASS::mvrnorm(n0, mu=true_params$a, Sigma=true_params$A)
X1 = MASS::mvrnorm(n1, mu=true_params$b, Sigma=true_params$B)
X = matrix(NA,nrow = n, ncol = d)
X[idx0,] = X0
X[idx1,] = X1
plot_cov(X,idx0,idx1,main=paste("d =",d))



# closed form
true_UOT = true_2Dnorm_UOT_map(true_params$a,true_params$b,
                               true_params$A,true_params$B, 
                               reg=reg, reg_marg=reg_m_kl)

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
a1 = rep(1/n1,n1)
a0 = rep(1/n0,n0)
ot <- import("ot")
maxiter = 2000
ubG = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl, div='kl',
                             numItermax=as.integer(maxiter))
ubG = ubG/sum(ubG)

# get empirical moments
get_emp_moments = function(ubG, X0, X1){
  d = ncol(X0)
  ubG_empE1 = rep(0, d*2)
  ubG_empE2 = matrix(0, 2*d, 2*d)
  t0 = Sys.time()
  for(i in 1:nrow(X0)){
    for(j in 1:nrow(X1)){
      ubG_empE1[1:d] = X0[i,]*ubG[j,i] + ubG_empE1[1:d]
      ubG_empE1[1:d+d] = X1[j,]*ubG[j,i] + ubG_empE1[1:d+d]
      
      ubG_empE2[1:d,1:d] = X0[i,]%*%t(X0[i,])*ubG[j,i] + ubG_empE2[1:d,1:d]
      ubG_empE2[1:d,1:d+d] = X0[i,]%*%t(X1[j,])*ubG[j,i] + ubG_empE2[1:d,1:d+d]
      ubG_empE2[1:d+d,1:d] = X1[j,]%*%t(X0[i,])*ubG[j,i] + ubG_empE2[1:d+d,1:d]
      ubG_empE2[1:d+d,1:d+d] = X1[j,]%*%t(X1[j,])*ubG[j,i] + ubG_empE2[1:d+d,1:d+d]
    }
  }
  t1 = Sys.time()
  elasped = t1 - t0
  return(list(ubG_empE1 = ubG_empE1, ubG_empE2 = ubG_empE2, elasped=elasped))
}

emp_ot_moments0 = get_emp_moments(ubG, X0, X1)

t1 = Sys.time()
elasped1 = t1 - t0

get_emp_moments_fast = function(ubG, X0, X1){
  d = ncol(X0)
  ubG_empE1 <- rep(0, d*2)
  ubG_empE2 <- matrix(0, 2*d, 2*d)
  t0 = Sys.time()
  # Calculate the empirical moments
  ubG_empE1[1:d] <- colSums(ubG %*% X0)
  ubG_empE1[(d+1):(2*d)] <- colSums(t(ubG) %*% X1)
  
  # Compute ubG_empE2
  for (j in 1:nrow(X1)) {
    weighted_X0 <- t(X0) * ubG[j, ]
    weighted_X1 <- t(X1[j, ]) * ubG[j, ]
    
    ubG_empE2[1:d, 1:d] <- ubG_empE2[1:d, 1:d] + weighted_X0 %*% X0
    ubG_empE2[1:d, (d+1):(2*d)] <- ubG_empE2[1:d, (d+1):(2*d)] + weighted_X0 %*% t(X1[j, ])
    ubG_empE2[(d+1):(2*d), 1:d] <- ubG_empE2[(d+1):(2*d), 1:d] + t(X1[j, ]) %*% weighted_X0
    ubG_empE2[(d+1):(2*d), (d+1):(2*d)] <- ubG_empE2[(d+1):(2*d), (d+1):(2*d)] + t(X1[j, ]) %*% t(weighted_X1)
  }
  t1 = Sys.time()
  elasped = t1 - t0
  return(list(ubG_empE1 = ubG_empE1, ubG_empE2 = ubG_empE2, elasped = elasped))
}

emp_ot_moments = get_emp_moments_fast(ubG, X0, X1)


