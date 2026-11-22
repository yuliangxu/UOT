# Single replication for callers supplying n_sim, d, true_params, reg and reg_m_kl.
figure1 <- new.env(parent = environment())
sys.source("simulation/ET_round2/MMD_heatmap.R", envir = figure1)

# Retain the original draw order for callers that set a seed before sourcing.
n0 <- n1 <- n_sim
n <- n0 + n1
idx0 <- sort(sample(seq_len(n), n0))
idx <- rep(1, n)
idx[idx0] <- 0
idx1 <- which(idx == 1)
X0 <- MASS::mvrnorm(n0, mu = true_params$a, Sigma = true_params$A)
X1 <- MASS::mvrnorm(n1, mu = true_params$b, Sigma = true_params$B)
X <- matrix(NA_real_, nrow = n, ncol = d)
X[idx0, ] <- X0
X[idx1, ] <- X1

population <- figure1$gaussian$gaussian_uot_population(
  true_params$a, true_params$b, true_params$A, true_params$B, reg, reg_m_kl
)
fit <- figure1$compute_uot_plan(X0, X1, reg, reg_m_kl)
ubG <- fit$G

start <- Sys.time()
mmd <- figure1$factorized_closed_form_mmd(X0, X1, ubG, population$mean, population$covariance)
out <- list(dist_cf = mmd$mmd2, elapsed_cf = as.numeric(difftime(Sys.time(), start, units = "secs")))
out
