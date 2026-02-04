# reproducible code for Figure 1 in ET submission


# One UOT matching example ------------------------------------------------

source("R/thresh_ot_func.R")
source("simulation/MR_functions.R")
library(reticulate)
use_python("~/Library/r-miniconda-arm64/bin/python", required = T) # choose your own python path
source_python("./python/sinkhorn_unbalanced_tv.py")

# curse of dim compared to other methods? --------------------------------------------

d = 2
n_sim = 10

# set reg params
reg_m_kl = 1e-2/5 # rho
reg = 1e-2 # epsilon



n_rep = 50
n_sim_grid = c(20, 40, 60, 80)
d_grid = c(2,3,4,5)

result_df = expand.grid(n_sim = n_sim_grid, rep = 1:n_rep, d = d_grid)
# result_df = expand.grid(n_sim = n_sim_grid, rep = 1:n_rep)
result_df$dist = rep(NA,nrow(result_df))
result_df$dist_cf = rep(NA,nrow(result_df))
result_df$time_emp = rep(NA,nrow(result_df))
result_df$time_cf = rep(NA,nrow(result_df))

ct = 0
for(d in d_grid){
  true_params = get_marg_norm(d)
  for(n_sim in n_sim_grid){
    print(paste("finished ",ct," out of ",nrow(result_df)," simulations"))
    for(rep in 1:n_rep){
      print(paste("n_sim:",n_sim,"rep:",rep,"d:",d))
      
      source("simulation/One_Rep_MMD.R")
      # result_df$dist[result_df$n_sim == n_sim & result_df$rep == rep & result_df$d == d] = out$dist
      result_df$dist_cf[result_df$n_sim == n_sim & result_df$rep == rep & result_df$d == d] = out$dist_cf
      # result_df$time_emp[result_df$n_sim == n_sim & result_df$rep == rep & result_df$d == d] = out$elapsed_emp
      result_df$time_cf[result_df$n_sim == n_sim & result_df$rep == rep & result_df$d == d] = out$elapsed_cf
      
      ct = ct+1
      saveRDS(result_df, paste("result/MR_R1Q2_1_sim",ct,"_total",nrow(result_df),".rds",sep=""))
      if(ct>1){
        unlink(paste("result/MR_R1Q2_1_sim",ct-1,"_total",nrow(result_df),".rds",sep=""))  
      }
      
    }
  }
  
}

head(result_df[result_df$d==2,])
head(result_df[result_df$d==4,])

result_df = readRDS(paste("result/MR_R1Q2_1_sim",800,"_total",800,".rds",sep=""))

# create a plot
library(ggplot2)
library(dplyr)
summary_df <- result_df %>%
  group_by(n_sim, d) %>%
  summarize(
    # mean_dist = median(dist),
    # lower_quantile = quantile(dist, 0.25),
    # upper_quantile = quantile(dist, 0.75),
    mean_dist_cf = median(dist_cf),
    lower_quantile_cf = quantile(dist_cf, 0.25),
    upper_quantile_cf = quantile(dist_cf, 0.75),
    # time_emp = mean(time_emp),
    time_cf = mean(time_cf)
  )
summary_df$d = as.factor(summary_df$d)
# ggplot(summary_df, aes(x = n_sim, y = mean_dist, color = d, group = d)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = lower_quantile, ymax = upper_quantile), width = 0.2) +
#   labs(
#     title = paste("Median and 25%,75% Quantiles of Empirical MMD over ", n_rep, " Replications",sep=""),
#     x = "Number of observations (N)",
#     y = "Distance (MMD)"
#   ) +
#   theme_minimal()

summary_df = summary_df[summary_df$d==2 | summary_df$d==3 | summary_df$d==5,]

MMD_plot = ggplot(summary_df, aes(x = sqrt(n_sim), y = sqrt(mean_dist_cf), color = d, group = d)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = sqrt(lower_quantile_cf), ymax = sqrt(upper_quantile_cf)), width = 0.2) +
  labs(
    title = paste("Median and 25%,75% Quantiles of closed-form MMD over ", n_rep, " Replications",sep=""),
    x = "Square-root Number of observations (sqrt(N))",
    y = "Distance (MMD)"
  ) +
  theme_minimal()

pdf("./plot/MMD_plot.pdf", width = 7, height = 4) # Adjust width and height as needed
print(MMD_plot)
dev.off()

# rate of convergence varying rho -----------------------------------------


