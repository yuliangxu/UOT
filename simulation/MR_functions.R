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


emp_MMD_fast = function(X,Y,X_weight){
  
  dist_matrix <- as.matrix(dist(X))
  kernel_matrix <- exp(-dist_matrix^2)
  weight_matrix <- outer(X_weight, X_weight)
  E_XX <- sum(weight_matrix * kernel_matrix)
  
  
  # Compute pairwise Euclidean distances between rows of X and Y
  dist_matrix <- as.matrix(dist(rbind(X, Y)))
  n <- nrow(X)
  m <- nrow(Y)
  dist_matrix <- dist_matrix[1:n, (n+1):(n+m)]
  
  # Compute the kernel values (Gaussian kernel in this case)
  kernel_matrix <- exp(-dist_matrix^2)
  
  # Compute the final value by summing the weighted kernel values
  E_XY <- sum(outer(X_weight, rep(1, m)) * kernel_matrix)
  
  
  dist_matrix <- as.matrix(dist(Y))
  kernel_matrix <- exp(-dist_matrix^2)
  E_YY <- sum(kernel_matrix)
  
  MMD = E_XX - 2*E_XY/nrow(Y) + E_YY/nrow(Y)^2
  
  
  return(MMD)
  
}

emp_MMD_fast2 = function(X,Y,X_weight){
  d = ncol(X)
  dist_matrix <- as.matrix(dist(X))
  kernel_matrix <- exp(-sqrt(dist_matrix^2/d))
  weight_matrix <- outer(X_weight, X_weight)
  E_XX <- sum(weight_matrix * kernel_matrix)
  
  
  # Compute pairwise Euclidean distances between rows of X and Y
  n <- nrow(X)
  m <- nrow(Y)
  dist_matrix_XY <- sqrt(outer(rowSums(X^2), rowSums(Y^2), `+`) - 2 * tcrossprod(X, Y))
  kernel_matrix_XY <- exp(-sqrt(dist_matrix_XY^2/d))
  E_XY <- sum(outer(X_weight, rep(1, m)) * kernel_matrix_XY)
  
  # Compute pairwise Euclidean distances for Y
  dist_matrix_Y <- as.matrix(dist(Y))
  kernel_matrix_Y <- exp(-sqrt(dist_matrix_Y^2/d))
  E_YY <- sum(kernel_matrix_Y)
  
  # Compute MMD
  MMD <- E_XX - 2 * E_XY / n + E_YY / (n * m)
  
  
  
  return(MMD)
  
}

closedform_MMD = function(X,X_weight,mu,Sigma){
  d = ncol(X)
  dist_matrix <- as.matrix(dist(X))
  kernel_matrix <- exp(-dist_matrix^2/d)
  weight_matrix <- outer(X_weight, X_weight)
  E_XX <- sum(weight_matrix * kernel_matrix)
  
  
  # E_XY
  joint_d = d
  H = 2*Sigma/joint_d + diag(rep(1,joint_d))
  H_inv = chol2inv(chol(H))
  H_det = det(H)
  quad_form = apply(X,1,function(x){
    t(x-mu) %*% H_inv %*% (x-mu)
  })
  E_XY = sum(exp(-quad_form/joint_d)*X_weight)/sqrt(H_det)
  
  
  
  # E_YY
  E_YY = 1.0/sqrt(det(Sigma/joint_d + 0.5*H))
  
  # Compute MMD
  MMD <- E_XX - 2 * E_XY  + E_YY 
  
  
  
  return(MMD)
  
}

ker = function(x,y){
  exp(-sum((x-y)^2))
}

generate_psd_matrix <- function(d,seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- .Random.seed  # Save the current seed
    set.seed(seed)            # Set the new seed
  }
  
  # Generate a random 8x8 matrix
  A <- matrix(rnorm(d*d), nrow=d, ncol=d)
  
  # Create a positive semi-definite matrix by multiplying A by its transpose
  PSD_matrix <- A %*% t(A)
  
  if (!is.null(seed)) {
    .Random.seed <<- old_seed  # Restore the previous seed
  }
  
  return(PSD_matrix)
}


get_marg_norm = function(d){
  d_char = as.character(d)
  out <- switch(d_char,
                "2" = list(
                    a = c(-0.4,0.1), # mean 0
                    b = c(0.5,-0.3), # mean 1
                    
                    A = matrix( c(0.7, -0.5,
                                  -0.5,  1), nrow=2, ncol=2), # var 0
                    
                    B = matrix( c(1, -0.34,
                                  -0.34,  0.7), nrow=2, ncol=2) # var 1)
                  ),
                "3" = {list(
                  a = c(-0.4,0.1,0.2), # mean 0
                  b = c(0.5,-0.3,-0.2), # mean 1
                  
                  A = generate_psd_matrix(3, seed = 2024),
                  
                  B = generate_psd_matrix(3, seed = 2024)
                )},
                "4" = {list(
                    a = c(-0.4,0.1,0.2,-0.3), # mean 0
                    b = c(0.5,-0.3,-0.2,0.3), # mean 1
                    
                    A = generate_psd_matrix(4, seed = 2024),
                    
                    B = generate_psd_matrix(4, seed = 2024)
                  )},
                "5" = {list(
                  a = c(-0.4,0.1,0.2,-0.3,0.1), # mean 0
                  b = c(0.5,-0.3,-0.2,0.3,-0.1), # mean 1
                  
                  A = generate_psd_matrix(5, seed = 2024),
                  
                  B = generate_psd_matrix(5, seed = 2024)
                )},
                "6" = {list(
                    a = c(-0.4,0.1,0.2,-0.3,0.1,-0.2), # mean 0
                    b = c(0.5,-0.3,-0.2,0.3,-0.1,0.2), # mean 1
                    
                    A = generate_psd_matrix(6, seed = 2024),
                    
                    B = generate_psd_matrix(6, seed = 2024)
                  )},
                "8" = {list(
                    a = c(-0.4,0.1,0.2,-0.3,0.1,-0.2,0.3,-0.1), # mean 0
                    b = c(0.5,-0.3,-0.2,0.3,-0.1,0.2,-0.3,0.1), # mean 1
                    
                    A = generate_psd_matrix(8, seed = 2024),
                    
                    B = generate_psd_matrix(8, seed = 2024)
                  )},
                "10" ={ list(
                    a = c(-0.4,0.1,0.2,-0.3,0.1,-0.2,0.3,-0.1, 0.7,-0.3), # mean 0
                    b = c(0.5,-0.3,-0.2,0.3,-0.1,0.2,-0.3,0.1, 0.3,-0.7), # mean 1
                    
                    A = generate_psd_matrix(10, seed = 2024),
                    
                    B = generate_psd_matrix(10, seed = 2024)
                  )}
  )
  
  return(out)
}

# Efficient computation of distance matrix
compute_distance_matrix <- function(X0, X1) {
  # Number of rows in each matrix
  n0 <- nrow(X0)
  n1 <- nrow(X1)
  
  # Compute the squared differences and sum them row-wise
  M <- sqrt(outer(rowSums(X0^2), rowSums(X1^2), `+`) - 2 * tcrossprod(X0, X1))
  
  return(M)
}