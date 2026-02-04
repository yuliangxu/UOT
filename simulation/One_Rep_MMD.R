## set input normal marginals ##
n0 = n_sim
n1 = n_sim
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
# plot_cov(X,idx0,idx1,main=paste("d =",d))



# closed form
true_UOT = true_2Dnorm_UOT_map(true_params$a,true_params$b,
                               true_params$A,true_params$B, 
                               reg=reg, reg_marg=reg_m_kl,d=d)

# distance matrix
# M = matrix(NA,n0,n1)
# for(i in 1:n0){
#   for(j in 1:n1){
#     M[i,j] = sqrt(sum((X0[i,] - X1[j,])^2))
#   }
# }

M = compute_distance_matrix(X0, X1)
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

# get MMD

X_cf = MASS::mvrnorm(n0*n1, mu = true_UOT$mu, Sigma = true_UOT$H)


# Generate indices for the combinations
i_indices <- rep(1:n0, each = n1)
j_indices <- rep(1:n1, times = n0)

# Fill X_joint and ubG_vec
X_joint <- cbind(X0[i_indices, ], X1[j_indices, ])
ubG_vec <- ubG[cbind(j_indices, i_indices)]


out = NULL
# t0 = Sys.time()
# out$dist = emp_MMD_fast2(X_joint, X_cf, ubG_vec)
# t1 = Sys.time()
# out$elapsed_emp = as.numeric(difftime(t1, t0),units = "secs")
# out

t0 = Sys.time()
out$dist_cf = closedform_MMD(X = X_joint, X_weight = ubG_vec, mu = true_UOT$mu, Sigma = true_UOT$H)
t1 = Sys.time()
out$elapsed_cf = as.numeric(difftime(t1, t0),units = "secs")
out
