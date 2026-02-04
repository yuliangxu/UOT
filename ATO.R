
library(reticulate)
source("./R/thresh_ot_func.R")
use_python("~/Library/r-miniconda-arm64/bin/python", required = T) # choose your own python path
source_python("./python/sinkhorn_unbalanced_tv.py")
# sim = 100
n0 = 100*sim
n1 = 10*sim
# d = 10
n = n1 + n0
n_rep = 100
# n_rep = 1

linear_bool = T
# covariate_case = 2

# what has worked:
# d = 6, sim = 50, case = 2, penalty = 1e-3/4
# d = 8, sim = 50, case = 2, penalty = 1e-3/4
# d = 10, sim = 50, case = 2, penalty = 1e-3/4

Y0_est = array(0,c(n1,3,n_rep)) # ot, knn truth
Y1_est = array(0,c(n0,3,n_rep))

output_all = NULL

for(iter in 1:n_rep){
  if(iter/10 - as.integer(iter/10) ==0){
    print(paste("iter = ", iter))
  } 
  
  
  idx0 = sort(sample(1:n,n0))
  idx = rep(1,n)
  idx[idx0] = 0
  idx1 = which(idx == 1)
  
  
  # par(mfrow = c(2,2))
  # case 1 (less overlap)
  if(covariate_case == 1){
    X = matrix(NA,nrow = n, ncol = d)
    X[idx0,] = rbind(matrix(rnorm(n0/2*d, mean = -1, sd = 2/sqrt(d)),ncol = d),
                     matrix(rnorm(n0/2*d, mean = -1, sd = 1/sqrt(d)),ncol = d))
    X[idx1,] = rbind(matrix(rnorm(n1/2*d, mean = 1, sd = 2/sqrt(d)),ncol = d),
                     matrix(rnorm(n1/2*d, mean = 1, sd = 1/sqrt(d)),ncol = d))
    # par(mfrow = c(1,2))
    # plot_cov(X,idx0,idx1,main="case 1 (less overlap)")
  }
  
  
  # case 2 (more overlap)
  if(covariate_case == 2){
    X = matrix(NA,nrow = n, ncol = d)
    X[idx0,] = rbind(matrix(rnorm(n0/2*d, mean = -1, sd = 0.5),ncol = d),
                     matrix(rnorm(n0/2*d, mean = 1/2, sd = 0.5),ncol = d))
    X[idx1,] = rbind(matrix(rnorm(n1/2*d, mean = 1/2, sd = 0.5),ncol = d),
                     matrix(rnorm(n1/2*d, mean = -1/2, sd = 0.5),ncol = d))
    # par(mfrow = c(1,1))
    # plot_cov(X,idx0,idx1,main="case 2 (more overlap)")
  }
  
  # case 3 (disc x conti: less overlap)
  if(covariate_case == 3){
    X = matrix(NA,nrow = n, ncol = d)
    X[idx0,] = rbind(matrix(rnorm(n0/2*d, mean = -1, sd = 2),ncol = d),
                     matrix(rnorm(n0/2*d, mean = 1/2, sd = 1),ncol = d))
    X[idx1,] = rbind(matrix(rnorm(n1/2*d, mean = 1, sd = 2),ncol = d),
                     matrix(rnorm(n1/2*d, mean = 1/2, sd = 1),ncol = d))
    X[,1] = as.integer(X[,1])
    # par(mfrow = c(1,1))
    # plot_cov(X,idx0,idx1,main="case 3 (disc x conti: less overlap)")
  }
  
  
  # case 4 (disc x conti: more overlap)
  if(covariate_case == 4){
    X = matrix(NA,nrow = n, ncol = d)
    X[idx0,] = rbind(matrix(rnorm(n0/2*d, mean = -1, sd = 0.5),ncol = d),
                     matrix(rnorm(n0/2*d, mean = 1/2, sd = 0.5),ncol = d))
    X[idx1,] = rbind(matrix(rnorm(n1/2*d, mean = 1/2, sd = 0.5),ncol = d),
                     matrix(rnorm(n1/2*d, mean = -1/2, sd = 0.5),ncol = d))
    X[,1] = as.integer(X[,1])
    # par(mfrow = c(1,1))
    # plot_cov(X,idx0,idx1,main="case 4 (disc x conti: overlap)")
  }
  
  
  
  
  Y = rep(NA,n)
  # 6d
  if(d == 6){
    if(linear_bool){
      Y0_all = -1 - apply(X[,1:3],1,sum) * apply(X[,4:6],1,sum) + 
        apply(X[,2:5],1,mean)+X[,4] + rnorm(n,sd = 1)
      Y1_all = 2 + 2*apply(X[,1:3],1,sum) + apply(X[,4:6],1,sum) +
        apply(X[,2:5],1,mean) - X[,4] + rnorm(n,sd = .5)
      # par(mfrow = c(2,3))
      # hist(Y1_all);hist(Y1_all[idx1]);hist(Y1_all[idx0]);
      # hist(Y0_all);hist(Y0_all[idx1]);hist(Y0_all[idx0]);
      # par(mfrow = c(1,1))
      # mean(Y1_all) - mean(Y0_all)
    }else{
      # nonlinear
      Y0_all = log((apply(X[,1:3],1,sum) +apply(X[,4:6],1,sum) )^2 + 0.5*exp(-X[,2]/10)) + rnorm(n,sd = 1)
      Y1_all = log(10+ exp(1-apply(X[,1:3],1,sum) - apply(X[,4:6],1,sum))) + rnorm(n,sd = .5)
    }
    
  }
  
  # 8d
  if(d == 8){
    if(linear_bool){
      Y0_all = -1 - apply(X[,1:4],1,sum) * apply(X[,5:8],1,sum) + 
        apply(X[,2:5],1,mean)+X[,4] + rnorm(n,sd = 1)
      Y1_all = 2 + 2*apply(X[,1:4],1,sum) + apply(X[,5:8],1,sum) +
        apply(X[,2:5],1,mean) - X[,4] + rnorm(n,sd = .5)
      # par(mfrow = c(2,3))
      # hist(Y1_all);hist(Y1_all[idx1]);hist(Y1_all[idx0]);
      # hist(Y0_all);hist(Y0_all[idx1]);hist(Y0_all[idx0]);
      # par(mfrow = c(1,1))
      # mean(Y1_all) - mean(Y0_all)
    }else{
      # nonlinear
      Y0_all = log((apply(X[,1:3],1,sum) +apply(X[,4:6],1,sum) )^2 + 0.5*exp(-X[,2]/10)) + rnorm(n,sd = 1)
      Y1_all = log(10+ exp(1-apply(X[,1:3],1,sum) - apply(X[,4:6],1,sum))) + rnorm(n,sd = .5)
    }
    
  }
  
  # 10d
  if(d == 10){
    if(linear_bool){
      Y0_all = -1 - apply(X[,1:5],1,sum) * apply(X[,6:10],1,sum) + apply(X[,1:5],1,sum)+X[,4] + rnorm(n,sd = 1)
      Y1_all = 2 + 2*apply(X[,1:4],1,sum) + apply(X[,5:8],1,sum) - X[,9]^2+ X[,10]^2 + rnorm(n,sd = .5)
    }else{
      # nonlinear
      Y0_all = log((apply(X[,1:3],1,sum) +apply(X[,4:8],1,sum) +X[,9]+X[,10])^2 + 0.5*exp(-X[,2]/10)) + rnorm(n,sd = 1)
      Y1_all = log(10+ exp(1-apply(X[,1:3],1,sum) - apply(X[,4:8],1,sum) -X[,9]-X[,10])) + rnorm(n,sd = .5)
      # par(mfrow = c(2,3))
      # hist(Y1_all);hist(Y1_all[idx1]);hist(Y1_all[idx0]);
      # hist(Y0_all);hist(Y0_all[idx1]);hist(Y0_all[idx0]);
      # par(mfrow = c(1,1))
      # mean(Y1_all) - mean(Y0_all)
    }
    
  }
  
  
  # 4d
  if(d == 4){
    if(linear_bool){
      Y0_all = -1 - X[,1] * X[,2] + X[,3]+X[,4] + rnorm(n,sd = 1)
      Y1_all = 2 + 2*X[,1] + X[,2] - X[,3]^2+ X[,4]^2 + rnorm(n,sd = .5)
    }else{
      # nonlinear
      Y0_all = log((X[,1]+X[,2]+X[,3]+X[,4])^2 + 0.5*exp(-X[,2]/10)) + rnorm(n,sd = 1)
      Y1_all = log(10+ exp(1-X[,1]- X[,2]-X[,3]-X[,4])) + rnorm(n,sd = .5)
      # par(mfrow = c(2,3))
      # hist(Y1_all);hist(Y1_all[idx1]);hist(Y1_all[idx0]);
      # hist(Y0_all);hist(Y0_all[idx1]);hist(Y0_all[idx0]);
      # par(mfrow = c(1,1))
      # mean(Y1_all) - mean(Y0_all)
    }
    
  }
  
  # 3d
  if(d == 3){
    if(linear_bool){
      Y0_all = -1 - X[,1] * X[,2] + X[,3] + rnorm(n,sd = 1)
      Y1_all = 2 + 2*X[,1] + X[,2] - X[,3]^2 + rnorm(n,sd = .5)
    }else{
      # nonlinear
      Y0_all = log((X[,1]+X[,2]+X[,3])^2 + 0.5*exp(-X[,2]/10)) + rnorm(n,sd = 1)
      Y1_all = log(10+ exp(1-X[,1]- X[,2]-X[,3])) + rnorm(n,sd = .5)
    }
    
  }
  
  
  # 2d
  if(d == 2){
    if(linear_bool){
      # linear 
      Y0_all = -1 - X[,1] * X[,2]  + rnorm(n,sd = 1)
      Y1_all = 2 + 2*X[,1] + X[,2] + rnorm(n,sd = .5)
    }else{
      # nonlinear
      Y0_all = log((X[,1]+X[,2])^2 + 0.5*exp(-X[,2]/10)) + rnorm(n,sd = 1)
      Y1_all = log(10+ exp(1-X[,1]- X[,2])) + rnorm(n,sd = .5)
    }
    
    
   
    # par(mfrow = c(2,3))
    # hist(Y1_all);hist(Y1_all[idx1]);hist(Y1_all[idx0]);
    # hist(Y0_all);hist(Y0_all[idx1]);hist(Y0_all[idx0]);
    # par(mfrow = c(1,1))
    # mean(Y1_all) - mean(Y0_all)
  }
  
  
  Y0_group = cbind(Y0_all[idx0], Y1_all[idx0])
  Y1_group = cbind(Y0_all[idx1], Y1_all[idx1])
  
  
  X0 = X[idx0,]
  X1 = X[idx1,]
  Y0 = Y0_all[idx0]
  Y1 = Y1_all[idx1]
  
  a1 = rep(1/n1,n1)
  a0 = rep(1/n0,n0)
  center0 = apply(X0, 2, mean)
  center1 = apply(X1, 2, mean)
  
  M = matrix(NA,n0,n1)
  for(i in 1:n0){
    for(j in 1:n1){
      M[i,j] = sqrt(sum((X0[i,] - X1[j,])^2))
    }
  }
  M = M/max(M)
  M = t(M)
  
  
  # plot_2hist(Y0_all,"Y0_all")
  
  # OT
  ot <- import("ot")
  if(d>4){
    # 1e-3/2 for LargeN3 (KNN1=TV), 1e-3/4 for largeN2(OT=TV)
    # 1e-3/3 for LargeN4 (KNN1=TV)
    # go with LargeN2 result: removes more outliers
    reg_m_kl = 1e-3/6
    reg = 1e-3/6
    maxiter = 2000
    G_ub_tv = sinkhorn_knopp_unbalanced_tv(a1, a0, M, reg, reg_m_kl,if_tv = 1,numItermax=as.integer(maxiter))
    summary(apply(G_ub_tv,1,sum))
    table(apply(G_ub_tv,1,sum)==0)
    ubG = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl, div='kl',numItermax=as.integer(maxiter))
    summary(apply(ubG,1,sum))
    table(apply(ubG,1,sum)==0)
  }else{
    reg_m_kl = 1e-3/5
    reg = 1e-3/5
    maxiter = 2000
  }
  
  ubG = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl, div='kl',numItermax=as.integer(maxiter))
  Y0_ot = rep(NA,n)
  Y0_ot[idx0] = Y0
  Y0_ot[idx1] = (ubG %*% Y0) / apply(ubG,1,sum)
  Y1_ot = rep(NA,n)
  Y1_ot[idx0] = (t(ubG) %*% Y1 ) / apply(ubG,2,sum)
  Y1_ot[idx1] = Y1
  
  
  
  # OT-tv
  G_ub_tv = sinkhorn_knopp_unbalanced_tv(a1, a0, M, reg, reg_m_kl,if_tv = 1,numItermax=as.integer(maxiter))
  # summary(apply(G_ub_tv,1,sum))
  Y0_tv = rep(NA,n)
  Y0_tv[idx0] = Y0
  Y0_tv[idx1] = (G_ub_tv %*% Y0) / apply(G_ub_tv,1,sum)
  Y1_tv = rep(NA,n)
  Y1_tv[idx0] = (t(G_ub_tv) %*% Y1 ) / apply(G_ub_tv,2,sum)
  # Y1_marg = apply(ubG,1,sum)/sum(ubG)
  Y1_tv[idx1] = Y1
  
  
  # OT-tvow
  if(d>4){
    Y0_tvow = Y0_tv;Y1_tvow = Y1_tv;
  }
  if(d %in% c(2,3,4)){
    pi0 = apply(G_ub_tv,2,sum)/sum(G_ub_tv)
    pi1 = apply(G_ub_tv,1,sum)/sum(G_ub_tv)
    Y0_tvow = Y0_tv;Y1_tvow = Y1_tv;
    if(any(pi0==0)){
      Y1_tvow[idx0[pi0==0]] = NA
    }
    if(any(pi1==1)){
      Y0_tvow[idx1[pi1==0]] = NA
    }
  }
  
  
  # OT2
  ot <- import("ot")
  reg_m_kl = 0.005
  reg = 0.005
  
  ubG2 = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl, div='kl')
  Y0_ot2 = rep(NA,n)
  Y0_ot2[idx0] = Y0
  Y0_ot2[idx1] = (ubG2 %*% Y0) / apply(ubG2,1,sum)
  Y1_ot2 = rep(NA,n)
  Y1_ot2[idx0] = (t(ubG2) %*% Y1 ) / apply(ubG2,2,sum)
  Y1_ot2[idx1] = Y1
  
  # OT3
  ot <- import("ot")
  reg_m_kl = 0.01
  reg = 0.01
  
  ubG3 = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl, div='kl')
  Y0_ot3 = rep(NA,n)
  Y0_ot3[idx0] = Y0
  Y0_ot3[idx1] = (ubG3 %*% Y0) / apply(ubG3,1,sum)
  Y1_ot3 = rep(NA,n)
  Y1_ot3[idx0] = (t(ubG3) %*% Y1 ) / apply(ubG3,2,sum)
  Y1_ot3[idx1] = Y1
  
  
  # KNN (k=3)
  knn3 = get_knn(M,3)
  Y0_knn3 = rep(NA,n); Y1_knn3 = rep(NA,n)
  Y0_knn3[idx0] = Y0; Y1_knn3[idx1] = Y1;
  Y0_knn3[idx1] = apply(knn3$X1_knn,1,function(x){mean(Y0[x])})
  Y1_knn3[idx0] = apply(knn3$X0_knn,1,function(x){mean(Y1[x])})
  
  # KNN (k=1)
  knn1 = get_knn(M,1)
  Y0_knn1 = rep(NA,n); Y1_knn1 = rep(NA,n)
  Y0_knn1[idx0] = Y0; Y1_knn1[idx1] = Y1;
  Y0_knn1[idx1] = apply(knn1$X1_knn,2,function(x){mean(Y0[x])})
  Y1_knn1[idx0] = apply(knn1$X0_knn,2,function(x){mean(Y1[x])})
  
  
  # ===== run ps weighting ======
  Tassign = rep(0,n)
  Tassign[idx1] = 1
  ps = glm( Tassign ~ X, family=binomial(link='logit'))
  ps_e =  cbind(1,X) %*% ps$coefficients
      # overlap weight (Fan Li)
      Y_ow = list(Y1 = Y1*(1-ps_e[idx1]),Y0 = Y0*(ps_e[idx0]))
      ATE_ow = mean(Y_ow$Y1, na.rm=T) - mean(Y_ow$Y0, na.rm=T)
  
  ps_e = 1/(1+exp(-ps_e))
  ps_e[ps_e<1e-2] = NA
  ps_e[ps_e>1-1e-2] = NA
  Y_ps = list(Y1 = Y1/ps_e[idx1],Y0 = Y0/(1-ps_e[idx0]))
  ATE_ps = mean(Y_ps$Y1, na.rm=T) - mean(Y_ps$Y0, na.rm=T)
  
  # overlap weight (Fan Li)
  Y_ow = list(Y1 = Y1*(1-ps_e[idx1]),Y0 = Y0*(ps_e[idx0]))
  ATE_ow = mean(Y_ow$Y1, na.rm=T) - mean(Y_ow$Y0, na.rm=T)
  
  # doubly robust
  
  if(d == 4){
    pred0 = lm(Y0 ~ X0[,1] + X0[,2] + X0[,3]+ X0[,4])
    pred1 = lm(Y1 ~ X1[,1] + X1[,2] + X1[,3]+ X1[,4])
  }
  
  if(d == 3){
    pred0 = lm(Y0 ~ X0[,1] + X0[,2] + X0[,3])
    pred1 = lm(Y1 ~ X1[,1] + X1[,2] + X1[,3])
  }
  
  if(d == 2){
    pred0 = lm(Y0 ~ X0[,1] + X0[,2])
    pred1 = lm(Y1 ~ X1[,1] + X1[,2])
    
  }
  
  
  
  if(d %in% c(2,3,4)){
    Y0_pred = cbind(1, X) %*% pred0$coefficients
    Y1_pred = cbind(1, X) %*% pred1$coefficients
    ATE_dr = mean(Y1_pred - Y0_pred) + mean(  (Y1-Y1_pred[idx1])/(ps_e[idx1]) - 
                                                (Y0-Y0_pred[idx0])/(1-ps_e[idx0]), na.rm = T )
    
  }else{
    ATE_dr = 0
  }
  
  
  # doubly robust oracle
  if(d == 4){
    pred0 = lm(Y0 ~ X0[,1]*X0[,2] + X0[,3]*X0[,4])
    pred1 = lm(Y1 ~ X1[,1]*X1[,2] + X1[,3]*X1[,4])
    Y0_pred = cbind(1, X, X[,1]*X[,2],X[,3]*X[,4]) %*% pred0$coefficients # change this for d=3
    Y1_pred = cbind(1, X, X[,1]*X[,2],X[,3]*X[,4]) %*% pred1$coefficients
  }
  if(d == 3){
    pred0 = lm(Y0 ~ X0[,1]*X0[,2] + X0[,3])
    pred1 = lm(Y1 ~ X1[,1]*X1[,2] + X1[,3])
    Y0_pred = cbind(1, X, X[,1]*X[,2]) %*% pred0$coefficients # change this for d=3
    Y1_pred = cbind(1, X,X[,1]*X[,2]) %*% pred1$coefficients
  }
  if(d == 2){
    pred0 = lm(Y0 ~ X0[,1]*X0[,2] )
    pred1 = lm(Y1 ~ X1[,1]+X1[,2] )
    Y0_pred = cbind(1, X, X[,1]*X[,2]) %*% pred0$coefficients
    Y1_pred = cbind(1, X[,1], X[,2]) %*% pred1$coefficients
  }
  
  if(d %in% c(2,3,4)){
    ATE_dr_oracle = mean(Y1_pred - Y0_pred) + mean(  (Y1-Y1_pred[idx1])/(ps_e[idx1]) - 
                                                       (Y0-Y0_pred[idx0])/(1-ps_e[idx0]), na.rm = T )
    
  }else{
    ATE_dr_oracle=0
  }
  
  
  
  ATE_OT = mean(Y1_ot - Y0_ot, na.rm = T)
  ATE_OT2 = mean(Y1_ot2 - Y0_ot2, na.rm = T)
  ATE_OT3 = mean(Y1_ot3 - Y0_ot3, na.rm = T)
  ATE_TV = mean(Y1_tv - Y0_tv, na.rm = T)
  ATE_TVow = mean(Y1_tvow - Y0_tvow, na.rm = T)
  ATE_knn3 = mean(Y1_knn3 - Y0_knn3, na.rm = T)
  ATE_knn1 = mean(Y1_knn1 - Y0_knn1, na.rm = T)
  ATE_true = mean(Y1_all - Y0_all, na.rm = T)
  output = c(ATE_TV = ATE_TV,
             ATE_TVow = ATE_TVow,
             ATE_OT = ATE_OT, 
             ATE_OT2 = ATE_OT2, 
             ATE_OT3 = ATE_OT3, 
             ATE_knn3 = ATE_knn3, 
             ATE_knn1 = ATE_knn1,
             ATE_ow = ATE_ow,
             ATE_dr = ATE_dr, 
             ATE_dr_oracle = ATE_dr_oracle,
             ATE_ps = ATE_ps, 
             ATE_true = ATE_true)
  
  output_all = rbind(output_all , output)
  
}

# saveRDS(list(ATE = output_all), paste("./sim/sim_all2_n1",n1, "_n0", n0,"_case",
#                                       covariate_case,"_d",d,".rds",sep=""))
# df = as.data.frame(cbind(output_all[1:11]-output_all[12], output_all[1:11]))
# rownames(df ) = label[1:11]
# df

n_method = dim(output_all)[2]-1
bias = abs(apply(output_all[,1:n_method]-output_all[,n_method+1], 2,mean))
MSE = apply(output_all[,1:n_method]-output_all[,n_method+1], 2,function(x){mean(x^2)})
VAR = apply(output_all[,1:n_method], 2,function(x){mean((x-mean(x))^2)})

# result = cbind(bias,VAR, MSE)
# result[order(MSE),]