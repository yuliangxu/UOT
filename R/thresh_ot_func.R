library(ggplot2)

#' @title true_2Dnorm_UOT_map
#' @description Compute true unbalanced optimal transport map for 2D normal marginals
#' @param a mean of normal distribution 0
#' @param b mean of normal distribution 1
#' @param A covariance matrix of normal distribution 0
#' @param B covariance matrix of normal distribution 1
#' @param reg regularization parameter for the joint penalty
#' @param reg_marg regularization parameter for the marginal penalty
true_2Dnorm_UOT_map = function(a,b,A,B, reg, reg_marg, d=2){
  sigma2 = reg/2
  gamma = reg_marg
  lambda = sigma2 + gamma/2
  tau = gamma/(2*sigma2+gamma)
  Id = diag(rep(1,d))
  X = A + B + lambda * Id
  X_inv = chol2inv(chol(X))
  A_tilde = gamma/2*(Id - lambda*chol2inv(chol(A+lambda*Id)))
  B_tilde = gamma/2*(Id - lambda*chol2inv(chol(B+lambda*Id)))
  C_temp = 1/tau * A_tilde %*% B_tilde + sigma2^2/4*Id
  C = expm::sqrtm(C_temp) - sigma2/2*Id
  
  mu = rbind(a + A %*% X_inv %*%(b-a), b - B %*% X_inv %*%(b-a))
  Id_C_lambda = Id + C/lambda
  H = rbind( cbind( Id_C_lambda %*% (A-A%*%X_inv%*%A), C + Id_C_lambda %*% A %*% X_inv %*% B),
             cbind( t(C) + t(Id_C_lambda)%*%B%*%X_inv%*%A, Id_C_lambda%*%( B - B%*% X_inv%*% B ) ) )
  m_pi = sigma2^(d/2*sigma2/gamma+sigma2)
  m_pi = m_pi * det(C) * sqrt(det(A_tilde%*%B_tilde)^tau/det(A%*%B))^(1/(tau+1))
  m_pi = m_pi * exp(-1/2*sum(-t(b-a)%*%X_inv%*%(b-a))/(tau+1))
  m_pi = m_pi * sqrt(det(C - 2/gamma*A_tilde%*%B_tilde))
  
  return(list(mu = mu, H = H, m_pi = m_pi))
}

#' @title title hist_to_dist
#' @description Compute empirical distribution from histogram
#' @param X observed sample points
hist_to_dist = function( X, plot=T, n_bins = 100){
  # Step 2: Define the bin boundaries for the 2D grid
  x_bins <- seq(min(X[,1]), max(X[,1]), length.out = n_bins) # Example bin boundaries for x
  y_bins <- seq(min(X[,2]), max(X[,2]), length.out = n_bins) # Example bin boundaries for y
  
  # Step 3: Calculate the center of each bin
  x_centers <- (x_bins[-length(x_bins)] + x_bins[-1]) / 2
  y_centers <- (y_bins[-length(y_bins)] + y_bins[-1]) / 2
  
  # Step 4: Count the number of sample points falling into each bin
  
  # x_hist = hist(X[,1], breaks = x_bins, plot = F)
  x_hist$counts
  counts <- table(
    cut(X[,1], breaks = x_bins, include.lowest = TRUE),
    cut(X[,2], breaks = y_bins, include.lowest = TRUE)
  )
  
  # Convert the table to a data frame
  counts_df <- as.data.frame(counts)
  colnames(counts_df) <- c("x_bin", "y_bin", "Count")
  
  # Convert bin levels to numeric centers
  counts_df$x_center <- x_centers[as.numeric(counts_df$x_bin)]
  counts_df$y_center <- y_centers[as.numeric(counts_df$y_bin)]
  
  # Step 5: Normalize these counts to obtain the empirical distribution
  total_samples <- sum(counts_df$Count)
  counts_df$Probability <- counts_df$Count / total_samples
  
  

  # Step 6: Visualize the empirical distribution using a heatmap
  if(plot){
    p <- ggplot(counts_df, aes(x = x_center, y = y_center, fill = Probability)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "blue") +
      labs(title = "Empirical Distribution Heatmap", x = "X Center", y = "Y Center") +
      theme_minimal()
    
    print(p)
  }
  
  
  
  return(counts_df)
}

UNormal_density = function(X, mu, Sigma, m_pi){
  d = ncol(X)
  Sigma_inv = chol2inv(chol(Sigma))
  
  single_density = function(X, mu, Sigma, Sigma_inv, m_pi){
    return(1/sqrt((2*pi)^d*det(Sigma)) * exp(-1/2*t(X-mu)%*%Sigma_inv%*%(X-mu)) * m_pi)
  }
  # return(1/sqrt((2*pi)^d*det(Sigma)) * exp(-1/2*rowSums((X-mu)%*%Sigma_inv*(X-mu))) * m_pi)
  return(apply(X, 1, function(x) single_density(x, mu, Sigma, Sigma_inv, m_pi)))
  
}


# # test
# X_test = X[1:2,]
# mu_test = c(1,2)
# Sigma_test = matrix(c(1,0.5,0.5,1),2,2)
# Sigma_inv_test = chol2inv(chol(Sigma_test))
# m_pi_test = 1
# UNormal_density(X_test, mu_test, Sigma_test, m_pi_test)
# 
# UNormal_density(t(X_test[1,]), mu_test, Sigma_test, m_pi_test)
# UNormal_density(t(X_test[2,]), mu_test, Sigma_test, m_pi_test)
# 
# X_in = X_test[2,]
# 1/sqrt((2*pi)^d*det(Sigma_test)) * exp(-1/2*t(X_in-mu_test)%*%Sigma_inv_test%*%(X_in-mu_test)) * m_pi_test

# Function to combine rows of two matrices
combine_rows <- function(mat1, mat2) {
  n1 <- nrow(mat1)
  n2 <- nrow(mat2)
  result <- matrix(0, nrow = n1 * n2, ncol = ncol(mat1) + ncol(mat2))
  
  row_idx <- 1
  for (i in 1:n1) {
    for (j in 1:n2) {
      result[row_idx, ] <- c(mat1[i, ], mat2[j, ])
      row_idx <- row_idx + 1
    }
  }
  
  return(result)
}


get_ot_knn = function(x,y,k=3){
  if(length(x)!=length(y)){
    print("Error: length does not match")
  }else{
    w_idx = order(x,decreasing = T)[1:k]
    return(sum(y[w_idx]*x[w_idx]/sum(x[w_idx])))
  }
}
count_matched = function(G){
  n0_matched = apply(G,2, function(x){sum(x>0)})
  n1_matched = apply(G,1, function(x){sum(x>0)})
  return(list(n0_matched = n0_matched, n1_matched = n1_matched))
}
plot_cov = function(x,idx0,idx1,dim_plot = 1:2,main="cov plot"){
  
  if(dim(x)[2] != 2){
    print(paste("plot dim",dim_plot))
    x = x[,dim_plot]
  }
  s = min(x); l = max(x)
    plot(x[idx0,1], x[idx0,2],col="blue",pch=1,cex=.8,xlab = "x[1]",
         ylab = "x[2]",xlim = c(s,l),ylim = c(s,l),main=main)
    points(x[idx1,1], x[idx1,2],col="red",pch=2,cex=.8)
    legend("topright",c("control","treated"),
           cex=.8,col=c("blue","red"), pch=c(1,2))
    
  
}
plot_cov_neighbor = function(x,idx0,idx1,idx_neighbor,dim_plot = 1:2,main="cov plot",cex = 1){
  
  if(dim(x)[2] != 2){
    print(paste("plot dim",dim_plot))
    x = x[,dim_plot]
  }
  s = min(x[idx_neighbor]); l = max(x[idx_neighbor])
  idx0 = intersect(idx0,idx_neighbor)
  idx1 = intersect(idx1,idx_neighbor)
  plot(x[idx0,1], x[idx0,2],col="blue",pch=1,cex=cex,xlab = "x[1]",
       ylab = "x[2]",xlim = c(s,l),ylim = c(s,l),main=main)
  points(x[idx1,1], x[idx1,2],col="red",pch=2,cex=cex)
  legend("topright",c("control","treated"),
         cex=1,col=c("blue","red"), pch=c(1,2))
  
  
}
plot_cov_w_measure = function(x,idx0,idx1,a0,a1,scale=5,dim_plot = 1:2){
 
  if(dim(x)[2] != 2){
    print(paste("plot dim",dim_plot))
    x = x[,dim_plot]
  }
  s = min(x); l = max(x)
    plot(x[idx0,1], x[idx0,2],col="blue",pch=1,cex= a0*scale,xlab = "x[1]",
         ylab = "x[2]",xlim = c(s,l),ylim = c(s,l))
    points(x[idx1,1], x[idx1,2],col="red",pch=2,cex= a1*scale)
    legend("topright",c("control","treated"),
           cex=.5,col=c("blue","red"), pch=c(1,2))
    
  
}
# plot_cov(X,idx0,idx1)
add_top_match = function(G,X,idx0,idx1,m_vec = 1,add_cov_plot = T,dim_plot=1:2,main="cov plot",
                         cex = 1,
                         idx_neighbor=NULL){
  if(add_cov_plot){
    if(is.null(idx_neighbor)){
      plot_cov(X,idx0,idx1,dim_plot=dim_plot,main)
    }else{
      plot_cov_neighbor(X,idx0,idx1,idx_neighbor,dim_plot=dim_plot,main,cex)
    }
    
  }
  top_G = sort(G,decreasing = T)[m_vec]
  top_G = top_G[top_G>0]
  idx = sapply(top_G, function(x,G){which(G == x, arr.ind = TRUE)}, G = G)
  print("weights of top nonzero pairs")
  print(top_G)
  top_pairs = matrix(NA,nrow = length(top_G),ncol=2)
  colnames(top_pairs) = c("idx1","idx0")
  for( i in 1:length(top_G)){
    points(X[idx1[idx[1,i]],1],X[idx1[idx[1,i]],2], pch = as.character(m_vec[i]),cex=cex)
    points(X[idx0[idx[2,i]],1], X[idx0[idx[2,i]],2], pch = as.character(m_vec[i]),cex=cex)
    pair = cbind(X[idx1[idx[1,i]],], X[idx0[idx[2,i]],])
    top_pairs[i,] = c( idx[1,i], idx[2,i])
    lines(pair[1,],pair[2,],cex=cex)
  }
  return(top_pairs)
}




add_top_match_w_measure = function(G,X,idx0,idx1,a0,a1,m_vec = 1,add_cov_plot = T,
                                   scale = 5){
  if(add_cov_plot){
    plot_cov_w_measure(X,idx0,idx1,a0,a1,scale = scale)
  }
  top_G = sort(G,decreasing = T)[m_vec]
  idx = sapply(top_G, function(x,G){which(G == x, arr.ind = TRUE)}, G = G)
  for( i in 1:length(m_vec)){
    points(X[idx1[idx[1,i]],1],X[idx1[idx[1,i]],2], pch = as.character(m_vec[i]))
    points(X[idx0[idx[2,i]],1], X[idx0[idx[2,i]],2], pch = as.character(m_vec[i]))
  }
  
}

plot_outcome = function(Y0_all, Y1_all ,idx0_in,idx1_in,
                        sort_by_assign = F,
                        plot_vert = T,
                        sub_idx = 0, main=NULL){
  s = min(c(Y0_all,Y1_all),na.rm = T)-0.05
  l = max(c(Y0_all,Y1_all),na.rm = T)+0.05
  x = 1:length(Y0_all)
  if(sub_idx[1] == 0){
    idx0 = idx0_in; idx1 = idx1_in
  }else{
    idx0 = intersect(idx0_in,sub_idx)
    idx1 = intersect(idx1_in,sub_idx)
  }
  if(sort_by_assign){
    Y1_all[1:length(idx1)] = Y1_all[idx1]
    Y0_all[1:length(idx1)] = Y0_all[idx1]
    Y1_all[1:length(idx0) + length(idx1)] = Y1_all[idx0]
    Y0_all[1:length(idx0) + length(idx1)] = Y0_all[idx0]
    idx1 = 1:length(idx1)
    idx0 = 1:length(idx0) + length(idx1)
    x = 1:(length(idx0) + length(idx1))
    sub_idx = 1:(length(idx0) + length(idx1))
  }
  plot(x[idx0], Y1_all[idx0], col="blue",pch=16,cex=.8,
       xlab = "individual", ylab = "outcome",
       xlim = c(min(x),max(x)),ylim = c(s,l),main = main)
  points(x[idx1], Y1_all[idx1],col="red",pch=17,cex=.8)
  points(x[idx0], Y0_all[idx0],col="blue",pch=1,cex=.8)
  points(x[idx1], Y0_all[idx1],col="red",pch=2,cex=.8)
  if(plot_vert){
    abline(v=x[sub_idx], col="grey",cex=.1)
  }
  
  legend("topright",c("Y(1),T=0","Y(1), T=1",
                      "Y(0), T=0","Y(0), T=1"),
         cex=.5,col=c("blue","red","blue","red"), pch=c(16,17,1,2))
  # print("Note: shape and color both represents T assignment, filled/not represents treatment effect.")
  
}
get_knn = function(M,k=3){
  X1_knn = t(apply(M, 1, function(x){order(x)[1:k]}))
  X0_knn = t(apply(M, 2, function(x){order(x)[1:k]}))
  return(list(X0_knn = X0_knn, X1_knn = X1_knn))
}

plot_bs = function(bs,ATE=1,breaks=10, qqplot=T,plot_sd = T){
  if(ATE==1){
    p1 <-  hist(bs$ATE[,1],breaks = breaks,plot = F)
    p2 <-  hist(bs$ATE[,2], breaks = breaks,plot = F)
    p3 <-  hist(bs$ATE[,3], breaks = breaks,plot = F)
    p4 <-  hist(bs$ATE[,4], breaks = breaks,plot = F)
    p5 <-  hist(bs$ATE[,5], breaks = breaks,plot = F)
    p6 <-  hist(bs$ATE[,6], breaks = breaks,plot = F)
    p7 <-  hist(bs$ATE[,7], breaks = breaks,plot = F)
  }else{
    p1 <-  hist(bs$ATT[,1], breaks = breaks,plot = F)
    p2 <-  hist(bs$ATT[,2], breaks = breaks,plot = F)
    p3 <-  hist(bs$ATT[,3], breaks = breaks,plot = F)
    p4 <-  hist(bs$ATT[,4], breaks = breaks,plot = F)
    p5 <-  hist(bs$ATT[,5], breaks = breaks,plot = F)
    p6 <-  hist(bs$ATT[,6], breaks = breaks,plot = F)
    p7 <-  hist(bs$ATT[,7], breaks = breaks,plot = F)
  }
  
  names_all = colnames(bs$ATT)
  par(mfrow = c(2,3))
  plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,10))  # first histogram
  plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,10), add=T)  # second
  legend("topright",legend = names_all[c(1,2)],cex=0.5,
         fill=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)))
  plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,10))  # first histogram
  plot( p3, col=rgb(1,1,0,1/4), xlim=c(0,10), add=T)  # second
  legend("topright",legend = names_all[c(1,3)],cex=0.5,
         fill=c(rgb(0,0,1,1/4), rgb(1,1,0,1/4)))
  plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,10))  # first histogram
  plot( p4, col=rgb(1,0,1,1/4), xlim=c(0,10), add=T)  # second
  legend("topright",legend = names_all[c(1,4)],cex=0.5,
         fill=c(rgb(0,0,1,1/4), rgb(1,0,1,1/4)))
  plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,10))  # first histogram
  plot( p5, col=rgb(1,0,1/2,1/4), xlim=c(0,10), add=T)  # second
  legend("topright",legend = names_all[c(1,5)],cex=0.5,
         fill=c(rgb(0,0,1,1/4), rgb(1,0,1/2,1/4)))
  plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,10))  # first histogram
  plot( p6, col=rgb(1,1/2,1,1/4), xlim=c(0,10), add=T)  # second
  legend("topright",legend = names_all[c(1,6)],cex=0.5,
         fill=c(rgb(0,0,1,1/4), rgb(1,1/2,1,1/4)))
  plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,10))  # first histogram
  plot( p7, col=rgb(1,1/3,1,1/4), xlim=c(0,10), add=T)  # second
  legend("topright",legend = names_all[c(1,7)],cex=0.5,
         fill=c(rgb(0,0,1,1/4), rgb(1,1/3,1,1/4)))
  par(mfrow = c(1,1))
  
  if(qqplot){
    par(mfrow = c(2,3))
    if(ATE==1){
      for(i in 1:6+1){
        plot(bs$ATE[,1], bs$ATE[,i],
             xlab = names_all[1], ylab = names_all[i],main = "ATE",asp=1)
        abline(0,1,col="red")
      }
    }else{
      for(i in 1:6+1){
        plot(bs$ATT[,1], bs$ATT[,i],
             xlab = names_all[1], ylab = names_all[i],main = "ATT",asp=1)
        abline(0,1,col="red")
      }
    }
  }
  if(plot_sd){
    for(i in 1:6+1){
      if(i %in% c(2,5,6)){
        plot(bs$ATT_sd[,1], bs$ATT_sd[,i],
             xlab = names_all[1], ylab = names_all[i],main = "ATT_sd",asp=1)
        abline(0,1,col="red")
      }else{
        plot.new()
      }
      
    }
  }
  par(mfrow = c(1,1))
}
plot_bs_unbalanced = function(bs,ATE=1,breaks=10, qqplot=T,plot_sd = T,hist_plot=F){
  names_all = colnames(bs$ATT)
  if(hist_plot){
    if(ATE==1){
      p1 <-  hist(bs$ATE[,1],breaks = breaks,plot = F)
      p2 <-  hist(bs$ATE[,2], breaks = breaks,plot = F)
      p3 <-  hist(bs$ATE[,3], breaks = breaks,plot = F)
      p4 <-  hist(bs$ATE[,4], breaks = breaks,plot = F)
      p5 <-  hist(bs$ATE[,5], breaks = breaks,plot = F)
      p6 <-  hist(bs$ATE[,6], breaks = breaks,plot = F)
    }else{
      p1 <-  hist(bs$ATT[,1], breaks = breaks,plot = F)
      p2 <-  hist(bs$ATT[,2], breaks = breaks,plot = F)
      p3 <-  hist(bs$ATT[,3], breaks = breaks,plot = F)
      p4 <-  hist(bs$ATT[,4], breaks = breaks,plot = F)
      p5 <-  hist(bs$ATT[,5], breaks = breaks,plot = F)
      p6 <-  hist(bs$ATT[,6], breaks = breaks,plot = F)
    }
    
    
    # par(mfrow = c(2,2))
    plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,10))  # first histogram
    plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,10), add=T)  # second
    legend("topright",legend = names_all[c(1,2)],cex=0.5,
           fill=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)))
    plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,10))  # first histogram
    plot( p3, col=rgb(1,1,0,1/4), xlim=c(0,10), add=T)  # second
    legend("topright",legend = names_all[c(1,3)],cex=0.5,
           fill=c(rgb(0,0,1,1/4), rgb(1,1,0,1/4)))
    plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,10))  # first histogram
    plot( p4, col=rgb(1,0,1,1/4), xlim=c(0,10), add=T)  # second
    legend("topright",legend = names_all[c(1,4)],cex=0.5,
           fill=c(rgb(0,0,1,1/4), rgb(1,0,1,1/4)))
    plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,10))  # first histogram
    plot( p5, col=rgb(1,0,1/2,1/4), xlim=c(0,10), add=T)  # second
    legend("topright",legend = names_all[c(1,5)],cex=0.5,
           fill=c(rgb(0,0,1,1/4), rgb(1,0,1/2,1/4)))
    plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,10))  # first histogram
    plot( p6, col=rgb(1,0,1/2,1/4), xlim=c(0,10), add=T)  # second
    legend("topright",legend = names_all[c(1,5)],cex=0.5,
           fill=c(rgb(0,0,1,1/4), rgb(1,0,1/2,1/4)))
  }
  
  # par(mfrow = c(1,1))
  
  if(qqplot){
    print("qqplot...")
    # par(mfrow = c(2,2))
    if(ATE==1){
      for(i in 1:5+1){
        plot(bs$ATE[,1], bs$ATE[,i],
             xlab = names_all[1], ylab = names_all[i],main = "ATE",asp=1)
        abline(0,1,col="red")
      }
    }else{
      for(i in 1:5+1){
        plot(bs$ATT[,1], bs$ATT[,i],
             xlab = names_all[1], ylab = names_all[i],main = "ATT",asp=1)
        abline(0,1,col="red")
      }
    }
  }
  if(plot_sd){
    print("plot_sd...")
    for(i in 1:5+1){
      if(i %in% c(2,4,5)){
        print(paste("i=",i))
        plot(bs$ATT_sd[,1], bs$ATT_sd[,i],
             xlab = names_all[1], ylab = names_all[i],main = "ATT_sd",asp=1)
        abline(0,1,col="red")
      }else{
        plot.new()
      }
      
    }
  }
  # par(mfrow = c(1,1))
}
