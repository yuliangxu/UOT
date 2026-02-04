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

## --- Precompute pieces you already have ---
X1 <- as.matrix(df_obs_std[df_obs_std[[Tvar]] == 1, covars])
X0 <- as.matrix(df_obs_std[df_obs_std[[Tvar]] == 0, covars])
n1 <- nrow(X1); n0 <- nrow(X0)
Y1 <- df_obs_std[df_obs_std[[Tvar]] == 1, Yvar]
Y0 <- df_obs_std[df_obs_std[[Tvar]] == 0, Yvar]

a1 <- rep(1/n1, n1)
a0 <- rep(1/n0, n0)

## Cost matrix (Euclidean). Faster & clearer than double for-loops.
D  <- as.matrix(dist(rbind(X1, X0)))
M  <- D[seq_len(n1), n1 + seq_len(n0)]
M  <- M / max(M)          # scale to [0,1] as you did

## --- Define your grids (edit these to your “given range”) ---
# rho_grid <- signif(seq(1e-4, 1e2, length.out = 20), 3)  # reg_m_kl
rho_grid <- signif(10^seq(-4, 1, length.out = 20), 3)# reg_m_kl
eps_grid <- signif(10^seq(-4, -2, length.out = 15), 3)  # reg



## Storage
res <- list()

## --- Grid search over (rho, eps) ---
k <- 0
for (rho in rho_grid) {
  for (eps in eps_grid) {
    k <- k + 1
    message(sprintf("#===== rho (reg_m_kl)=%.3g, eps (reg)=%.3g", rho, eps))
    
    ubG3 <- ot$sinkhorn_unbalanced(a1, a0, M, eps, rho, div='kl')
    uot_dist_compare = evaluate_pairwise_G(X1, X0, ubG3, treat_weighting = "by_rowmass")
    
    
    
    att_G <- att_from_G(Y1, Y0, ubG3, X0,
                        G_thresh = 0,
                        row_weighting = "rowmass",return_matched_covariates = TRUE)  # matches the formula
    
    
    X0_matched <- att_G$X0_matched
    X0_matched_original <- as.data.frame(X0_matched)
    names(X0_matched_original) <- covars
    ## Undo standardization for continuous vars (as in your code)
    X0_matched_original[cont_vars] <- sweep(
      X0_matched_original[cont_vars], 2, scaling_params$sds, "*"
    )
    X0_matched_original[cont_vars] <- sweep(
      X0_matched_original[cont_vars], 2, scaling_params$means, "+"
    )
    X0_matched_original[true_cont_vars] <- round(X0_matched_original[true_cont_vars])
    
    summary_tbl <- hellinger_summary(
      nsw_t = nsw_t,
      X0_matched_original = X0_matched_original,
      nsw_c = nsw_c,
      covars = covars,
      cont_vars = true_cont_vars
    )
    
    H_sum_val <- colSums(summary_tbl[, c("H_matched", "H_exp")])[1]
    
    matched_controls_per_treat <- apply(ubG3, 1, function(x) sum(x > 0))
    med_matches <- median(matched_controls_per_treat)
    
    res[[k]] <- data.frame(
      rho = rho,
      eps = eps,
      H_sum = H_sum_val,
      median_matched_controls = med_matches,
      UOT_att = att_G$att
    )
  }
}

## Final summary table
summary_table <- do.call(rbind, res)
row.names(summary_table) <- NULL

## Optional: sort by a criterion, e.g., smallest H_sum
summary_table <- summary_table[order(summary_table$H_sum), ]

summary_table
# -------- 5) visualize (median_matched_controls)
library(dplyr)
library(ggplot2)
summary_table$median_matched_controls = factor(summary_table$median_matched_controls)
ggplot(summary_table,
       aes(x = factor(rho),
           y = factor(eps),
           fill = median_matched_controls)) +
  geom_tile(color = "white") +
  scale_fill_viridis_d(
    option = "C",
    name = "Median matched controls"
  ) +
  labs(
    x = expression(rho~"(reg_m_kl)"),
    y = expression(epsilon~"(reg)"),
    title = "Median matched controls across (rho, eps)"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# -------- 5) visualize (UOT_att)
library(dplyr)
library(ggplot2)

target <- 886.3
best_row <- summary_table |>
  filter(is.finite(UOT_att)) |>
  slice_min(abs(UOT_att - target), n = 1)

ggplot(summary_table, aes(x = factor(rho), y = factor(eps), fill = UOT_att)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(
    option = "C",
    # limits = c(-3000, 3000),   # ← clip here
    # oob = scales::squish,
    name = "UOT_att"
  ) +
  # highlight: draw a border around the best tile
  geom_tile(
    data = best_row,
    aes(x = factor(rho), y = factor(eps)),
    fill = NA, color = "red", linewidth = 1.3
  ) +
  # optional label on top of the highlighted tile
  geom_text(
    data = best_row,
    aes(label = paste0("≈", round(UOT_att, 1))),
    color = "red", 
    size=3,
    # fontface = "bold",
    vjust = -0.1
  )+
  labs(x = expression(rho~"(reg_m_kl)"),
       y = expression(epsilon~"(reg)"),
       fill = "UOT_att",
       title = "UOT_att across (rho, eps)\nHighlighted closets to 886.3") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- plot compare to ATT_exp
library(dplyr)
library(ggplot2)

target <- ATT_exp

df <- summary_table %>%
  filter(is.finite(UOT_att)) %>%
  mutate(
    rho_f = factor(rho, levels = sort(unique(rho))),
    eps_f = factor(eps, levels = sort(unique(eps))),
    err   = abs(UOT_att - target)
  )

min_err <- min(df$err, na.rm = TRUE)
tol <- 1e-10   # relax if needed, e.g. 1e-8

best_rows <- df %>%
  filter(abs(err - min_err) <= tol)

ggplot(df, aes(x = rho_f, y = eps_f, fill = UOT_att)) +
  geom_tile() +
  scale_fill_gradient2(
    low  = "darkblue",
    mid  = "white",
    high = "darkred",
    midpoint = 0,
    name = "UOT_att"
  ) +
  # highlight exactly one tile
  geom_tile(
    data = best_rows,
    aes(x = rho_f, y = eps_f),
    fill = NA, color = "red", linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = best_rows,
    aes(x = rho_f, y = eps_f, label = sprintf("≈%.1f", UOT_att)),
    color = "red", fontface = "bold", size = 3,
    inherit.aes = FALSE
  ) +
  labs(
    title = "UOT_att across (rho, eps)\nHighlighted closest to ATT_exp",
    x = expression(rho~"(reg_m_kl)"),
    y = expression(epsilon~"(reg)")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# ---- inspect matrix version of the heatmap
library(dplyr)
library(tidyr)

UOT_mat_df <- summary_table %>%
  filter(is.finite(UOT_att)) %>%
  mutate(
    rho_f = factor(rho, levels = sort(unique(rho))),
    eps_f = factor(eps, levels = sort(unique(eps)))
  ) %>%
  select(eps_f, rho_f, UOT_att) %>%
  pivot_wider(names_from = rho_f, values_from = UOT_att) %>%
  arrange(eps_f)

View(UOT_mat_df)
