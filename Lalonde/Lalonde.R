## ========= Setup =========
need <- c("cobalt", "sandwich", "lmtest", "utils")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, quiet = TRUE)

library(cobalt)
library(sandwich)
library(lmtest)

source("./Lalonde/lalonde_help.R")

base <- "https://users.nber.org/~rdehejia/data/"



## ========= 1) Load randomized NSW (Lalonde) =========
nsw_t  <- read.table(url(paste0(base, "nsw_treated.txt")), header = FALSE)
names(nsw_t) <- c("treat","age","educ","black","hisp","married","nodegree","re75","re78")

nsw_c  <- read.table(url(paste0(base, "nsw_control.txt")), header = FALSE)
names(nsw_c) <- c("treat","age","educ","black","hisp","married","nodegree","re75","re78")

# Sanity checks
stopifnot(all(nsw_t$treat == 1), all(nsw_c$treat == 0))

## ========= 2) Experimental benchmark ATT =========
Yvar <- "re78"
Tvar <- "treat"
ATT_exp <- mean(nsw_t[[Yvar]]) - mean(nsw_c[[Yvar]])

## ========= 3) Build observational sample (NSW treated + PSID controls) =========
psid1_c <- read_nsw(paste0(base, "psid_controls.txt"), has_re74 = TRUE)
psid2_c <- read_nsw(paste0(base, "psid2_controls.txt"), has_re74 = TRUE)
psid3_c <- read_nsw(paste0(base, "psid3_controls.txt"), has_re74 = TRUE)
psid_c <- rbind(psid1_c, psid2_c, psid3_c)

stopifnot(all(psid_c$treat == 0))
names(psid_c) = c("treat","age","educ","black","hisp","married","nodegree","re74","re75","re78")
psid_c$re74 = NULL

# df_obs <- rbind(nsw_t, nsw_c, psid_c)  # treated from NSW, controls from PSID
df_obs <- rbind(nsw_t, psid_c)  # treated from NSW, controls from PSID
covars <- c("age","educ","black","hisp","married","nodegree","re75")  # predictors used in DW
# cont_vars <- c("age", "educ", "re75")
cont_vars <- covars
true_cont_vars <- c("age", "educ", "re75")

df_obs_std <- df_obs
means <- colMeans(df_obs_std[cont_vars])
sds   <- apply(df_obs_std[cont_vars], 2, sd)
scaling_params <- list(means = means, sds = sds)
# Replace with standardized values
df_obs_std[cont_vars] <- sweep(df_obs_std[cont_vars], 2, scaling_params$means, "-")
df_obs_std[cont_vars] <- sweep(df_obs_std[cont_vars], 2, scaling_params$sds, "/")


X_std = df_obs_std[cont_vars]
X_original <- sweep(X_std, 2, scaling_params$sds, "*")
X_original <- sweep(X_original, 2, scaling_params$means, "+")

X_conti_cov = df_obs[cont_vars]
print(paste0("recon err=", mean(rowSums((as.matrix(X_original) - as.matrix(X_conti_cov))^2)))) 

## ========= 4) UOT =========

library(reticulate)
source("./R/thresh_ot_func.R")
use_python("~/Library/r-miniconda-arm64/bin/python", required = T) # choose your own python path
source_python("./python/sinkhorn_unbalanced_tv.py")
# UOT with reg = 0.01
ot <- import("ot")

# reg_m_kl = 1e-3
# reg = 1e-3
X1 = as.matrix(df_obs_std[df_obs_std[[Tvar]] == 1,covars])
X0 = as.matrix(df_obs_std[df_obs_std[[Tvar]] == 0,covars])
n1 = nrow(X1)
n0 = nrow(X0)
Y1 = df_obs_std[df_obs_std[[Tvar]] == 1,Yvar]
Y0 = df_obs_std[df_obs_std[[Tvar]] == 0,Yvar]
n = n1+n0

a1 = rep(1/n1,n1)
a0 = rep(1/n0,n0)

M = matrix(NA,n0,n1)
for(i in 1:n0){
  for(j in 1:n1){
    M[i,j] = sqrt(sum((X0[i,] - X1[j,])^2))
  }
}
M = M/max(M)
M = t(M)

reg_m_kl = 1e-3
reg = 5e-3

H_sum = NULL
UOT_att_list = NULL
median_matched_controls = NULL
rho_list = NULL
eps_list = NULL
reg_lsit = seq(1e-3,1e-4,length.out = 10)
for(reg in reg_lsit){
  reg_m_kl = reg
  reg = 5e-4
  print(paste("#================ reg=",reg))
  ubG3 = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl, div='kl')
  rho_list = c(rho_list, reg_m_kl)
  eps_list = c(eps_list, reg)
  idx0 = which(df_obs_std[[Tvar]] == 0)
  idx1 = which(df_obs_std[[Tvar]] == 1)
  
  
  
  uot_dist_compare = evaluate_pairwise_G(X1, X0, ubG3, treat_weighting = "by_rowmass")
  
  
  
  att_G <- att_from_G(Y1, Y0, ubG3, X0,
                      G_thresh = 0,
                      row_weighting = "rowmass",return_matched_covariates = TRUE)  # matches the formula
  
  print(paste0("UOT:","ATT = ",round(att_G$att)))
  UOT_att_list = c(UOT_att_list, att_G$att)
  print(paste0("Exp benchmark ATT = ",round(ATT_exp)))
  
  
  X0_matched = att_G$X0_matched
  X0_matched_original = as.data.frame(X0_matched)
  names(X0_matched_original) = covars
  X0_matched_original[cont_vars] <- sweep(X0_matched_original[cont_vars] , 2, scaling_params$sds, "*")
  X0_matched_original[cont_vars]  <- sweep(X0_matched_original[cont_vars] , 2, scaling_params$means, "+")
  X0_matched_original[true_cont_vars] = round(X0_matched_original[true_cont_vars])
  
  
  round(colMeans(X0_matched_original),digit=2)
  round(colMeans(nsw_t[covars]),digit=2)
  round(colMeans(nsw_c[covars]),digit=2)
  
  # covars
  # chosen_cov = covars[1] # age
  # chosen_cov = covars[7] # age
  # 
  # par(mfrow = c(1,length(cont_vars)))
  # for(chosen_cov in cont_vars){
  #   plot(nsw_t[[chosen_cov]],  X0_matched_original[[chosen_cov]],asp=1, 
  #        xlab = "treatment", ylab = "UOT matched control",
  #        main = paste0("UOT matched covariate: ",chosen_cov)
  #   );abline(0,1)
  # }
  
  discrete_vars = setdiff(covars, cont_vars)
  
  summary_tbl <- hellinger_summary(
    nsw_t = nsw_t,
    X0_matched_original = X0_matched_original,
    nsw_c = nsw_c,
    covars = covars,
    cont_vars = true_cont_vars
  )
  
  # print(summary_tbl, digits = 3)
  H_sum = c(H_sum,colSums(summary_tbl[,c("H_matched","H_exp")])[1])
  matched_controls_per_treat = apply(ubG3, 1, function(x){sum(x>0)})
  median_matched_controls = c(median_matched_controls, median(matched_controls_per_treat))
  results = cbind(rho = rho_list, eps = eps_list ,H_sum=H_sum, median_matched_controls = median_matched_controls,UOT_att = UOT_att_list)
  
}

cbind(rho = rho_list, eps = eps_list ,H_sum=H_sum, median_matched_controls = median_matched_controls,UOT_att = UOT_att_list)


overlap_hist(nsw_t[[chosen_cov]], nsw_c[[chosen_cov]],  X0_matched_original[[chosen_cov]],
             labels = c(paste0(chosen_cov,":Treatment1"), "NSW control", "UOT_matched control"),
             normalize = "density",
             add_density = TRUE, bins=40)



# check number of match pairs
reg_m_kl = 1
reg =  6e-4
thresh = 0
ubG3 = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl, div='kl')
matched_controls_per_treat = apply(ubG3, 1, function(x){sum(x>thresh)})
summary(matched_controls_per_treat)
print(paste0("range of G:(",min(ubG3),",",max(ubG3),")"))
r_k = rowSums(ubG3)
Y0_mapped = (ubG3%*%Y0) /r_k
diff = Y1-Y0_mapped
weighted_diff = sum(diff*r_k) / sum(r_k);weighted_diff
par(mfrow=c(1,2))
r_k_norm = r_k/sum(r_k)
thresh = 1e-4
plot(diff[r_k_norm>thresh],r_k_norm[r_k_norm>thresh],main=thresh)


# --------- UOT with fixed eps and rho ------------
reg = 6e-4
reg_m_kl = 6e-4
ubG3 = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl, div='kl')
att_G <- att_from_G(Y1, Y0, ubG3, X0,
                    G_thresh = 0,
                    row_weighting = "rowmass",return_matched_covariates = TRUE)  # matches the formula
uot_dist_compare = evaluate_pairwise_G(X1, X0, ubG3, treat_weighting = "by_rowmass")
UOT_att = att_G$att

# otG = ot$sinkhorn(a1, a0, M, reg, numItermax=as.integer(1000)) # balanced
# matched_controls_per_treat = apply(otG, 1, function(x){sum(x>thresh)})
# summary(matched_controls_per_treat)
# print(paste0("range of G:(",min(otG),",",max(otG),")"))


# ----------------- Competing method: KNN ------------------
# KNN (k=3)
knn3 = get_knn(M,3)
Y0_knn3 = rep(NA,n); Y1_knn3 = rep(NA,n)
Y0_knn3[idx0] = Y0; Y1_knn3[idx1] = Y1;
Y0_knn3[idx1] = apply(knn3$X1_knn,1,function(x){mean(Y0[x])})
Y1_knn3[idx0] = apply(knn3$X0_knn,1,function(x){mean(Y1[x])})

ATT_knn3 = mean(Y1_knn3[idx1] - Y0_knn3[idx1])
ATE_knn3 = mean(Y1_knn3 - Y0_knn3)
print(paste0("KNN (k=3):","ATT = ",ATT_knn3,"; ATE = ",ATE_knn3))

# ------------- compare distribution ------------------- #

library(MatchIt)
# library(ebalance)
library(CBPS)
library(sbw)
library(cem)
library(WeightIt)

# unadjusted
ATT_unadjusted = with(df_obs, mean(re78[treat == 1]) - mean(re78[treat == 0]) )

# convenience: function for weighted mean difference (ATT)
weighted_mean <- function(y, w) sum(y * w) / sum(w)
ATT_diff <- function(y, t, w) {
  mean(y[t == 1]) - weighted_mean(y[t == 0], w[t == 0])
}

# Propensity score nearest neighbor
m_nearest <- matchit(treat ~ age + educ + black + hisp + married + nodegree + re75,
                     caliper = 0.1,
                     data = df_obs, method = "nearest")

d_nearest <- match.data(m_nearest)
ATT_nearest <- with(d_nearest, mean(re78[treat == 1]) - mean(re78[treat == 0])); ATT_nearest

# Covariate Balancing Propensity Score (CBPS)
cbps_fit <- CBPS(treat ~ age + educ + black + hisp + married + nodegree + re75,
                 data = df_obs,ATT  = TRUE)
w_cbps <- cbps_fit$weights
ATT_cbps <- ATT_diff(df_obs$re78, df_obs$treat, w_cbps); ATT_cbps


# Mahalanobis matching (covariate distance)
m_maha <- matchit(treat ~ age + educ + black + hisp + married + nodegree + re75,
                  data = df_obs, method = "nearest", distance = "mahalanobis")
d_maha <- match.data(m_maha)
ATT_maha <- with(d_maha, mean(re78[treat == 1]) - mean(re78[treat == 0])); ATT_maha

# Stable Balancing Weights (SBW)
res_sbw <- sbw_att(dat = df_obs, ind = "treat", out = "re78", covars = covars,
                   bal_std = "target",
                   bal_tol = 0.02, bal_alg = FALSE)

res_sbw$att  


# Overlap weights (WeightIt)
w_glm <- weightit(treat ~ age + educ + black + hisp + married + nodegree + re75,
                  data = df_obs, method = "glm", estimand = "ATE")  # logistic PS
ps <- w_glm$ps
w_overlap <- with(df_obs, ps * (1 - ps))
w_overlap[df_obs$treat == 1] <- w_overlap[df_obs$treat == 1] / sum(w_overlap[df_obs$treat == 1])
w_overlap[df_obs$treat == 0] <- w_overlap[df_obs$treat == 0] / sum(w_overlap[df_obs$treat == 0])

ATT_overlap <-ATT_diff(df_obs$re78, df_obs$treat, w_overlap)
ATT_overlap


# summarize in a table
method_names = c("Unadjusted","PS Nearest","CBPS","Mahalanobis","SBW","Overlap Weights","UOT","Experimental")
ATT_values = c(ATT_unadjusted, ATT_nearest, ATT_cbps, ATT_maha, res_sbw$att, ATT_overlap,UOT_att,ATT_exp)
results_table = data.frame(Method = method_names, ATT = ATT_values)
print(results_table)
knitr::kable(results_table, format = "latex")


aggregate(. ~ df_obs$treat, data = df_obs[,true_cont_vars], summary)
discrete_vars = setdiff(covars, true_cont_vars)
aggregate(. ~ df_obs$treat, data = df_obs[,discrete_vars], function(x){mean(x)})


