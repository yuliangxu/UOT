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
reg_m_kl = 1e-2
reg = 1e-2
maxiter = 2000
ubG = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl, div='kl',
                             numItermax=as.integer(maxiter))


# compare with the truth
true_UOT = true_2Dnorm_UOT_map(a,b,A,B, reg=reg, reg_marg=reg_m_kl)
true_UOT$m_pi

true_UOT$H

jointX = combine_rows(X0, X1) # expand X0 row 1 and all rows in X1
true_UOT$density = UNormal_density(jointX,true_UOT$mu, true_UOT$H, true_UOT$m_pi)

true_UOT_df = data.frame(x=jointX[,1], y=jointX[,2], density = true_UOT$density/sum(true_UOT$density))

ggplot(true_UOT_df, aes(x = x, y = y, fill = density)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "UOT density", x = "X Center", y = "Y Center") +
  theme_minimal()

true_UOT_df$sh_UOT = c(ubG)

ggplot(true_UOT_df, aes(x = x, y = y, fill = sh_UOT)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "UOT density", x = "X Center", y = "Y Center") +
  theme_minimal()


plot(true_UOT_df$sh_UOT, 
     true_UOT$density, 
     asp=1,
     xlab = "Sinkhorn UOT", ylab = "True UOT", main = "Sinkhorn UOT vs True UOT");abline(0,1)

# map between hist and dist -----------------------------------------------

n_side = 20
grids_lim = range(apply(X,2,range))
grids = BayesGPfit::GP.generate.grids(d=d,num_grids = n_side,grids_lim = grids_lim)

X1_dist = hist_to_dist(X1)
X0_dist = hist_to_dist(X0)

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
reg_m_kl = 1e-2
reg = 1e-2
maxiter = 2000
ubG = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl, div='kl',numItermax=as.integer(maxiter))


# sort matrix -------------------------------------------------------------

# Create a sample matrix A
A <- matrix(c(4, 2, 7, 1, 5, 6, 3, 8, 9), nrow = 3, byrow = TRUE)
print("Original matrix:")
print(A)

# Flatten the matrix into a vector
A_vector <- as.vector(A)

# Get the order of the sorted vector
sorted_order <- order(A_vector)

# Sort the vector
sorted_A_vector <- A_vector[sorted_order]

# Get the row and column indices of the sorted elements
row_indices <- ((sorted_order - 1) %% nrow(A)) + 1
col_indices <- ((sorted_order - 1) %/% nrow(A)) + 1

# Recreate the sorted matrix
sorted_matrix <- matrix(sorted_A_vector, nrow = nrow(A), ncol = ncol(A), byrow = FALSE)

# Display the results
print("Sorted matrix:")
print(sorted_matrix)
print("Row indices of sorted elements:")
print(row_indices)
print("Column indices of sorted elements:")
print(col_indices)


# subsampling -------------------------------------------------------------

# row 1-2 is X0, row 3-4 is X1
X_sinkhorn = matrix(nrow = n_sample, ncol = d*2)
ubG_norm1 = ubG/sum(ubG)
A_vector <- as.vector(ubG_norm1)

# Get the order of the sorted vector
sorted_order <- order(A_vector)

# Sort the vector
sorted_A_vector <- A_vector[sorted_order]

# Get the row and column indices of the sorted elements
row_indices <- ((sorted_order - 1) %% nrow(A)) + 1
col_indices <- ((sorted_order - 1) %/% nrow(A)) + 1

# Recreate the sorted matrix
sorted_matrix <- matrix(sorted_A_vector, nrow = nrow(A), ncol = ncol(A), byrow = FALSE)

# create n_sample
probs = runif(n_sample)

