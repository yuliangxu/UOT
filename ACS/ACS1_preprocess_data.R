source("./ACS/ACS_help.R")
batch_path = "./data"
seed_i = 1 # change this for each batch

# # use that env
# # Sys.setenv(RETICULATE_PYTHON = "/hpc/home/yx306/.virtualenvs/r-venv/bin/python")
# # use_virtualenv("r-venv", required = TRUE)
# library(reticulate)
# 
# # reticulate::use_python("/usr/bin/python3", required = TRUE)
# # reticulate::py_config()
# source("./R/thresh_ot_func.R")
# # use_python("~/Library/r-miniconda-arm64/bin/python", required = T) # choose your own python path
# # virtualenv_create("r-venv", python = "/usr/bin/python3")
# # virtualenv_install("r-venv", packages = c("POT"))
# 
library(reticulate)
source("./R/thresh_ot_func.R")
use_python("~/Library/r-miniconda-arm64/bin/python", required = T) # choose your own python path
source_python("./python/sinkhorn_unbalanced_tv.py")
# UOT with reg = 0.01
ot <- import("ot")


source_python("./python/sinkhorn_unbalanced_tv.py")

# 1) read and preprocess data --------------------------------------------------
df <- read.csv(gzfile("./data/usa_00002.csv.gz"))

couple_df <- make_spouse_pairs(df,
                         covars = c("AGE","EDUC","INCWAGE","OCC","IND","EMPSTAT","MARST"),
                         id_vars = c("YEAR"))

male_df <- couple_df$male_cov
female_df <- couple_df$female_cov
male_df <- male_df %>% rename_with(~ sub("^M_", "", .x))
female_df <- female_df %>% rename_with(~ sub("^F_", "", .x))

dim(male_df); dim(female_df)   # same number of rows (couples), numeric matrices for OT
head(male_df); head(female_df) 
# remove YEAR: all 2024 data
male_df$YEAR = NULL; female_df$YEAR = NULL
# remove suspicious data
bad_id <- "37178_4_5"
male_df   <- male_df   %>% dplyr::filter(couple_id != bad_id)
female_df <- female_df %>% dplyr::filter(couple_id != bad_id)
rm(df, couple_df) # to save memory

# ---- choose one covariate for 1D OT / binning ----
var1d <- "AGE"   # or "INCWAGE"
K <- 30          # number of bins

# pull global vectors
x1_all <- male_df[[var1d]]
x0_all <- female_df[[var1d]]

# optional: log-transform income globally (recommended)
if (var1d == "INCWAGE") {
  x1_all <- log1p(pmax(x1_all, 0))
  x0_all <- log1p(pmax(x0_all, 0))
}

# global breaks from pooled support
rng <- range(c(x1_all, x0_all), finite = TRUE)
breaks <- seq(rng[1], rng[2], length.out = K + 1)

# bin centers (useful later for barycentric map)
bin_centers <- 0.5 * (breaks[-1] + breaks[-length(breaks)])


# 2) create bootstrap samples (one batch)---------------------------------------------
batch <- bootstrap_couples_one(male_df, female_df, n_pairs = 1000, 
                               seed = seed_i)
cont_vars = c("AGE","INCWAGE")
# treat all covariates as continuous
male_proc   <- prep_ot_mat(batch$male,   
                           cont_vars = cont_vars)
female_proc <- prep_ot_mat(batch$female, 
                           cont_vars = cont_vars)

# 3) UOT (one batch)------------------------------------------------------------------

# ---- bootstrap sample: extract 1D covariate ----
x1 <- batch$male[,var1d]
x0 <- batch$female[,var1d]

if (var1d == "INCWAGE") {
  x1 <- log1p(pmax(x1, 0))
  x0 <- log1p(pmax(x0, 0))
}


# assign to global bins
bin1 <- cut(x1, breaks = breaks, include.lowest = TRUE, labels = FALSE)
bin0 <- cut(x0, breaks = breaks, include.lowest = TRUE, labels = FALSE)

# weighted empirical marginals on bins (length K)
a1_bin <- rep(0, K)
a0_bin <- rep(0, K)

# aggregate weights by bin (fast + NA-safe)
tmp1 <- tapply(batch$w_male, bin1, sum)
tmp0 <- tapply(batch$w_female, bin0, sum)

a1_bin[as.integer(names(tmp1))] <- as.numeric(tmp1)
a0_bin[as.integer(names(tmp0))] <- as.numeric(tmp0)

# normalize again to get empirical distributions on bins (recommended)
a1_bin <- a1_bin / sum(a1_bin)
a0_bin <- a0_bin / sum(a0_bin)

# cost matrix between bin centers (1D)
bin_centers <- 0.5 * (breaks[-1] + breaks[-length(breaks)])
M_bin <- abs(outer(bin_centers, bin_centers, "-"))
M_bin <- M_bin / max(M_bin)

reg_m_kl <- 10
reg      <- 1e-2

G_bin <- ot$sinkhorn_unbalanced(a1_bin, a0_bin, M_bin, reg, reg_m_kl, div = "kl")

# View(G_bin)
library(ggplot2)

K <- nrow(M_bin)
df_hm <- expand.grid(i = 1:K, j = 1:K)
df_hm$val <- as.vector(G_bin)

ggplot(df_hm, aes(x = j, y = i, fill = val)) +
  geom_tile() +
  scale_y_reverse() +
  labs(
    x = "Target bin (female)",
    y = "Source bin (male)",
    fill = "value",
    title = sprintf("Heatmap (eps = %.3g, rho = %.3g)", reg, reg_m_kl)
  ) +
  theme_minimal()

