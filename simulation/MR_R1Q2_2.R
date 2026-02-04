# One UOT matching example ------------------------------------------------

source("R/thresh_ot_func.R")
source("simulation/MR_functions.R")
library(reticulate)
use_python("~/Library/r-miniconda-arm64/bin/python", required = T) # choose your own python path
source_python("./python/sinkhorn_unbalanced_tv.py")

# curse of dim compared to other methods? --------------------------------------------

d = 2
true_params = get_marg_norm(d)
n_sim = 100

# set reg params
reg_m_kl = 1e-2/5 # rho
reg = 1e-2 # epsilon



n_rep = 50
rho_grid = c(0.001, 0.002, 0.0005, 0.01)

n0 = n_sim
n1 = n_sim
n = n0+n1
idx0 = sort(sample(1:n,n0))
idx = rep(1,n)
idx[idx0] = 0
idx1 = which(idx == 1)

result_df = expand.grid(rep = 1:n_rep, rho = c(Inf, rho_grid))
result_df$num_within_radius = rep(NA,nrow(result_df))
result_df$num_match_within_radius = rep(NA,nrow(result_df))
result_df$proportion = rep(NA,nrow(result_df))

ct = 0
for(i_rep in 1:n_rep){
  X0 = MASS::mvrnorm(n0, mu=true_params$a, Sigma=true_params$A)
  X1 = MASS::mvrnorm(n1, mu=true_params$b, Sigma=true_params$B)
  
  
  
  
  # closed form
  true_UOT = true_2Dnorm_UOT_map(true_params$a,true_params$b,
                                 true_params$A,true_params$B, 
                                 reg=reg, reg_marg=reg_m_kl,d=d)
  
  M = compute_distance_matrix(X0, X1)
  M = M/max(M)
  M = t(M)
  a1 = rep(1/n1,n1)
  a0 = rep(1/n0,n0)
  ot <- import("ot")
  maxiter = 5000
  
  # balanced OT
  otG = ot$sinkhorn(a1, a0, M, reg, numItermax=as.integer(maxiter))
  
  out_g = count_match_within_radius(otG,M, radius = 0.1, measure_lower=1e-2)
  result_df$num_within_radius[result_df$rep == i_rep & result_df$rho == Inf] = out_g[1]
  result_df$num_match_within_radius[result_df$rep == i_rep & result_df$rho == Inf] = out_g[2]
  result_df$proportion[result_df$rep == i_rep & result_df$rho == Inf] = out_g[3]
  
  print(paste("simulation ",ct,"finished out of ",n_rep*length(rho_grid)))
  for(reg_m_kl in rho_grid){
    
    print(paste("rho=",reg_m_kl,", i_rep = ",i_rep))
    # UOT estimator with uniform measure
    ubG = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl, div='kl',
                                 numItermax=as.integer(maxiter))
    out_ub = count_match_within_radius(ubG, M, radius = 0.1, measure_lower=1e-2)
    
    result_df$num_within_radius[result_df$rep == i_rep & result_df$rho == reg_m_kl] = out_ub[1]
    result_df$num_match_within_radius[result_df$rep == i_rep & result_df$rho == reg_m_kl] = out_ub[2]
    result_df$proportion[result_df$rep == i_rep & result_df$rho == reg_m_kl] = out_ub[3]
    
    ct = ct+1
  }
}

head(result_df)

# analyze -----------------------------------------------------------------


# create a plot
library(ggplot2)
library(dplyr)
summary_df <- result_df %>%
  group_by(rho) %>%
  summarize(
    mean_prop = median(proportion, na.rm = TRUE),
    lower = quantile(proportion, 0.25, na.rm = TRUE),
    upper = quantile(proportion, 0.75, na.rm = TRUE)
  )

summary_df$rho = as.factor(summary_df$rho)

summary_df

ggplot(summary_df, aes(x = rho, y = mean_prop)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(
    title = "Proportion of matched pairs within radius 0.1, with 25% and 75% quantiles bars, over 50 replications",
    x = "rho",
    y = "Proportion of matched pairs"
  ) +
  theme_minimal()
