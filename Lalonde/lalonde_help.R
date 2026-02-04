library(ggplot2)

# Helper to read whitespace-delimited text with known columns
read_nsw <- function(url_txt,has_re74 = TRUE) {
  df <- read.table(url(url_txt), header = FALSE)
  if (has_re74) {
    names(df) <- c("treat","age","educ","black","hisp","married","nodegree","re74","re75","re78")
  } else {
    names(df) <- c("treat","age","educ","black","hisp","married","nodegree","re75","re78")
  }
  df
}

UOT_impute = function(Y0,Y1,ubG3,idx0,idx1){
  # Y(0)
  Y0_ot3 = rep(NA,n)
  Y0_ot3[idx0] = Y0
  Y0_ot3[idx1] = (ubG3 %*% Y0) / apply(ubG3,1,sum)
  # Y(1)
  Y1_ot3 = rep(NA,n)
  Y1_ot3[idx0] = (t(ubG3) %*% Y1 ) / apply(ubG3,2,sum)
  Y1_ot3[idx1] = Y1

  # ATT_UOT = mean(Y1_ot3[idx1] - Y0_ot3[idx1])
  # ATE_UOT = mean(Y1_ot3 - Y0_ot3)
  return(list(Y1_ot3=Y1_ot3,Y0_ot3=Y0_ot3))
}


att_from_G <- function(Y1, Y0, G,
                       G_thresh = 0,
                       row_weighting = c("rowmass","equal"),
                       X0 = NULL,              # optional: n0 × p control covariates
                       return_matched_covariates = TRUE) {
  row_weighting <- match.arg(row_weighting)
  G = pmax(G, G_thresh)  # thresholding small values
  G  <- as.matrix(G); Y1 <- as.numeric(Y1); Y0 <- as.numeric(Y0)
  n1 <- length(Y1);   n0 <- length(Y0)
  stopifnot(nrow(G) == n1, ncol(G) == n0)
  if (any(G < 0, na.rm = TRUE)) stop("G must be nonnegative")
  
  r    <- rowSums(G)           # row masses r_i
  keep <- r > 0                # treated rows actually matched
  if (!any(keep)) stop("All treated rows have zero mass in G.")
  
  Gk  <- G[keep, , drop = FALSE]
  rk  <- r[keep]
  Y1k <- Y1[keep]
  
  ## ----- ATT under G -----
  # Per-treated synthetic control outcome: sum_j (G_ij / r_i) * Y0_j
  Y0_syn_k <- as.numeric((Gk %*% Y0) / rk)
  
  # Per-treated differences
  d <- Y1k - Y0_syn_k
  
  if (row_weighting == "rowmass") {
    # ATT_G = (1 / sum_i r_i) * sum_i r_i * (Y1i - sum_j G_ij / r_i * Y0j)
    att <- sum(rk * d) / sum(rk)
  } else {
    att <- mean(d)
  }
  
  ## ----- Matched covariates X0_matched (rowmass scheme) -----
  X0_matched <- NULL
  if (!is.null(X0) && return_matched_covariates) {
    X0m <- matrix(NA_real_, n1, ncol(X0))     # preserve original treated order
    colnames(X0m) <- colnames(X0)
    
    # Row-normalize G on kept rows, then multiply by X0
    W <- Gk
    W <- W / rk                                  # broadcast row-wise
    X0_syn_cov <- W %*% as.matrix(X0)            # n1_kept × p synthetic covariates
    
    X0m[keep, ] <- X0_syn_cov
    X0_matched  <- X0m
  }
  
  out <- list(
    att   = att,
    keep  = keep,            # which treated rows had positive mass
    row_mass = r             # r_i for reference
  )
  if (!is.null(X0_matched)) out$X0_matched <- X0_matched
  return(out)
}

# Inputs:
#   X1: n1 x p matrix/data.frame of covariates for treated
#   X0: n0 x p matrix/data.frame of covariates for controls
#   G : n1 x n0 pairwise weights (nonnegative; rows/cols may sum to zero)
# Options:
#   treat_weighting = "equal" (default) or "by_rowmass"
# Returns: list with SMDs, KS, per-treated fit, ESS, and drop fractions
evaluate_pairwise_G <- function(X1, X0, G, treat_weighting = c("equal","by_rowmass")) {
  treat_weighting <- match.arg(treat_weighting)
  X1 <- as.matrix(X1); X0 <- as.matrix(X0); G <- as.matrix(G)
  n1 <- nrow(X1); n0 <- nrow(X0); p <- ncol(X1)
  
  stopifnot(ncol(X1) == ncol(X0), nrow(G) == n1, ncol(G) == n0)
  if (any(G < 0, na.rm = TRUE)) stop("G must be nonnegative")
  
  # Keep only units with positive mass in G
  rmass <- rowSums(G)  # mass assigned to each treated
  cmass <- colSums(G)  # mass assigned to each control
  keep1 <- rmass > 0
  keep0 <- cmass > 0
  
  X1k <- X1[keep1, , drop = FALSE]
  X0k <- X0[keep0, , drop = FALSE]
  Gk  <- G[keep1, keep0, drop = FALSE]
  rmk <- rowSums(Gk)  # row masses in kept set
  cmk <- colSums(Gk)
  
  # ---- 1) Per-treated normalization (row-normalized W) ----
  # Synthetic control for treated i: \tilde X0_i = sum_j W_ij X0_j
  W <- Gk
  nz <- rmk > 0
  W[nz, ] <- W[nz, , drop = FALSE] / rmk[nz]
  # Per-treated synthetic covariates (n1_kept x p)
  X0_syn <- W %*% X0k
  
  # How to average across treated?
  if (treat_weighting == "equal") {
    w_treat <- rep(1, nrow(X1k))
  } else { # by_rowmass
    w_treat <- rmk
  }
  w_treat <- w_treat / sum(w_treat)
  
  mean_treat_row <- as.numeric(colSums(X1k * w_treat))
  mean_ctrl_row  <- as.numeric(colSums(X0_syn * w_treat))
  
  sd_ref <- apply(X1k, 2, sd)  # treated SD as reference (common in matching)
  sd_ref[sd_ref == 0] <- 1
  SMD_row <- (mean_treat_row - mean_ctrl_row) / sd_ref
  
  # KS distance per covariate using per-treated synthetic data
  KS_row <- sapply(seq_len(p), function(k) {
    # Compare empirical distribution of treated X1k[,k] vs synthetic X0_syn[,k]
    suppressWarnings(ks.test(X1k[,k], X0_syn[,k])$statistic)
  })
  
  # Per-treated fit (Euclidean) — how close each treated's synthetic control is
  per_unit_fit <- sqrt(rowSums((X1k - X0_syn)^2))
  fit_summary  <- c(mean = mean(per_unit_fit), median = median(per_unit_fit),
                    q90 = quantile(per_unit_fit, 0.9), max = max(per_unit_fit))
  
  # ---- 2) Global aggregation (column-sum weights on controls) ----
  w_ctrl_global <- cmk
  w_ctrl_global <- w_ctrl_global / sum(w_ctrl_global)
  mean_ctrl_global <- as.numeric(colSums(X0k * w_ctrl_global))
  
  # Treated mean to compare against (same treated subset; choose equal weighting)
  mean_treat_global <- colMeans(X1k)
  SMD_global <- (mean_treat_global - mean_ctrl_global) / sd_ref
  
  # ---- Weight concentration diagnostics ----
  ESS_global <- (sum(Gk)^2) / sum(Gk^2)    # overall ESS of the transport plan
  # Entropy per row: dispersion of weights used to construct each treated's synthetic control
  row_entropy <- -rowSums(W * log(pmax(W, 1e-12)))
  entropy_summary <- c(mean = mean(row_entropy), median = median(row_entropy))
  
  # ---- Dropped mass / units ----
  frac_treated_dropped <- mean(!keep1)
  frac_control_dropped <- mean(!keep0)
  mass_treated_dropped <- sum(rmass[!keep1]) / sum(rmass + (rmass==0))  # robust denom
  mass_control_dropped <- sum(cmass[!keep0]) / sum(cmass + (cmass==0))
  
  list(
    # Main balance outputs
    SMD_row      = SMD_row,
    KS_row       = KS_row,
    SMD_global   = SMD_global,
    # Per-treated synthetic fit
    per_unit_fit = per_unit_fit,
    fit_summary  = fit_summary,
    # Concentration diagnostics
    ESS_global   = ESS_global,
    row_entropy  = row_entropy,
    entropy_summary = entropy_summary,
    # Dropped info
    frac_treated_dropped = frac_treated_dropped,
    frac_control_dropped = frac_control_dropped,
    mass_treated_dropped = mass_treated_dropped,
    mass_control_dropped = mass_control_dropped,
    # Means reported (handy for tables)
    mean_treat_row      = mean_treat_row,
    mean_ctrl_row       = mean_ctrl_row,
    mean_treat_global   = mean_treat_global,
    mean_ctrl_global    = mean_ctrl_global
  )
}
overlap_hist <- function(...,
                         labels = NULL,
                         bins = NULL, binwidth = NULL,
                         normalize = c("count", "density", "probability"),
                         colors = NULL,
                         alpha = 0.45,
                         add_density = TRUE,
                         na_rm = TRUE) {
  normalize <- match.arg(normalize)
  samples <- list(...)
  stopifnot(length(samples) >= 2)
  # all numeric
  if (!all(vapply(samples, is.numeric, TRUE)))
    stop("All samples must be numeric.")
  
  # drop NAs if requested
  if (na_rm) samples <- lapply(samples, function(v) v[is.finite(v)])
  
  n <- length(samples)
  if (is.null(labels)) labels <- paste0("Sample", seq_len(n))
  
  df <- data.frame(
    value = unlist(samples, use.names = FALSE),
    group = factor(rep(labels, lengths(samples)), levels = labels)
  )
  
  rng <- range(df$value, na.rm = TRUE)
  
  # Shared binning
  if (is.null(bins) && is.null(binwidth)) {
    fd_bw <- function(v) 2 * (IQR(v, na.rm = TRUE) / (length(v)^(1/3)))
    bw <- max(vapply(samples, fd_bw, numeric(1)))
    if (!is.finite(bw) || bw <= 0) bw <- diff(rng) / 30
    binwidth <- bw
  }
  
  # Colors
  if (is.null(colors)) {
    if (n <= 3) {
      colors <- c("#1f77b4", "#ff7f0e", "#2ca02c")[seq_len(n)]
    } else {
      colors <- scales::hue_pal()(n)
    }
  }
  
  suppressPackageStartupMessages(require(ggplot2))
  
  p <- ggplot(df, aes(x = value, fill = group))
  
  if (normalize == "count") {
    p <- p + geom_histogram(aes(y = after_stat(count)),
                            position = "identity",
                            bins = bins, binwidth = binwidth,
                            alpha = alpha, color = "white")
    y_lab <- "Count"
  } else if (normalize == "density") {
    p <- p + geom_histogram(aes(y = after_stat(density)),
                            position = "identity",
                            bins = bins, binwidth = binwidth,
                            alpha = alpha, color = "white")
    y_lab <- "Density"
  } else {
    p <- p + geom_histogram(aes(y = after_stat(count / sum(count))),
                            position = "identity",
                            bins = bins, binwidth = binwidth,
                            alpha = alpha, color = "white")
    y_lab <- "Probability"
  }
  
  p <- p +
    scale_fill_manual(values = colors) +
    coord_cartesian(xlim = rng, expand = FALSE) +
    labs(x = NULL, y = y_lab, fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")
  
  if (add_density && normalize != "count") {
    p <- p +
      geom_density(aes(x = value, y = after_stat(density), color = group),
                   linewidth = 0.8, inherit.aes = FALSE, data = df) +
      scale_color_manual(values = colors) +
      guides(color = "none")
  }
  
  p
}

w_normalize <- function(w) {
  if (is.null(w)) return(NULL)
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  s <- sum(w)
  if (s == 0) rep(0, length(w)) else w / s
}

wquantile <- function(x, w, probs) {
  # Weighted type-1 quantile (step function)
  stopifnot(length(x) == length(w))
  o <- order(x); x <- x[o]; w <- w[o]
  w <- w / sum(w)
  cw <- cumsum(w)
  sapply(probs, function(p) x[which(cw >= p)[1]])
}
wasserstein1_1d <- function(x, y, wx = NULL, wy = NULL, grid = 1024L) {
  keepx <- is.finite(x) & !is.na(x); wx <- if (is.null(wx)) NULL else wx[keepx]
  keepy <- is.finite(y) & !is.na(y); wy <- if (is.null(wy)) NULL else wy[keepy]
  if (length(x) == 0L || length(y) == 0L) return(NA_real_)
  wx <- w_normalize(if (is.null(wx)) rep(1, length(x)) else wx)
  wy <- w_normalize(if (is.null(wy)) rep(1, length(y)) else wy)
  u <- (seq_len(grid) - 0.5) / grid
  qx <- wquantile(x, wx, u)
  qy <- wquantile(y, wy, u)
  mean(abs(qx - qy))  # integral |F^{-1}_x - F^{-1}_y|
}
# Make a weighted probability vector over a specified union of levels
prob_from_labels <- function(labels, w = NULL, levels_union = NULL) {
  labs <- as.character(labels)
  if (is.null(levels_union)) levels_union <- sort(unique(labs))
  if (is.null(w)) w <- rep(1, length(labs)) else w <- pmax(0, as.numeric(w))
  
  # normalize weights within the observed labels
  s <- sum(w[is.finite(w)])
  if (s == 0) w[] <- 0 else w <- w / s
  
  # tabulate onto the full union (fill missing with 0)
  p <- setNames(numeric(length(levels_union)), levels_union)
  if (length(labs)) {
    tt <- tapply(w, labs, sum)
    p[names(tt)] <- as.numeric(tt)
  }
  # already sums to 1 by construction
  p
}

# KL helper (nats); returns 0 if a cell is zero in both
KL <- function(p, q) {
  idx <- p > 0 & q > 0
  if (!any(idx)) return(0)
  sum(p[idx] * log(p[idx] / q[idx]))
}

# Categorical distance with explicit alignment
cat_distance <- function(a, b, wa = NULL, wb = NULL,
                         metric = c("tv", "hellinger", "js", "sqrtjs"),
                         levels_union = NULL) {
  metric <- match.arg(metric)
  a_chr <- as.character(a); b_chr <- as.character(b)
  
  # Build the shared level union and align both p, q to it (same length & order)
  if (is.null(levels_union))
    levels_union <- sort(unique(c(a_chr, b_chr)))
  
  pa <- prob_from_labels(a_chr, wa, levels_union)
  pb <- prob_from_labels(b_chr, wb, levels_union)
  
  if (metric == "tv") {
    # total variation in [0,1]
    return(0.5 * sum(abs(pa - pb)))
  }
  if (metric == "hellinger") {
    # Hellinger distance in [0,1]
    return(sqrt(0.5 * sum((sqrt(pa) - sqrt(pb))^2)))
  }
  # JS and sqrt(JS)
  m <- 0.5 * (pa + pb)
  js <- 0.5 * KL(pa, m) + 0.5 * KL(pb, m)
  if (metric == "js") return(js)
  sqrt(js)
}
# Hellinger distance via shared histogram bins
hellinger_hist <- function(x, y, wx = NULL, wy = NULL,
                           breaks = "FD",
                           include_lowest = TRUE, right = TRUE,
                           return_squared = FALSE) {
  # keep finite
  x <- x[is.finite(x) & !is.na(x)]
  y <- y[is.finite(y) & !is.na(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  
  # shared breaks from pooled data (or user-supplied numeric vector)
  all_vals <- c(x, y)
  br <- if (is.character(breaks) && length(breaks) == 1) {
    hist(all_vals, breaks = breaks, plot = FALSE)$breaks
  } else {
    as.numeric(breaks)
  }
  if (length(br) < 2) return(NA_real_)
  
  # normalize weights to sum 1 within each sample
  nrmw <- function(w, n) {
    if (is.null(w)) rep(1/n, n) else {
      w <- pmax(0, as.numeric(w))
      s <- sum(w); if (s == 0) rep(1/n, n) else w / s
    }
  }
  wx <- nrmw(wx, length(x))
  wy <- nrmw(wy, length(y))
  
  # bin to probabilities on the SAME breaks
  bx <- cut(x, br, include.lowest = include_lowest, right = right)
  by <- cut(y, br, include.lowest = include_lowest, right = right)
  px <- as.numeric(tapply(wx, bx, sum)); if (anyNA(px)) px[is.na(px)] <- 0
  py <- as.numeric(tapply(wy, by, sum)); if (anyNA(py)) py[is.na(py)] <- 0
  
  # ensure they sum to 1 (guards against rare numeric drift)
  if ((sx <- sum(px)) > 0) px <- px / sx
  if ((sy <- sum(py)) > 0) py <- py / sy
  
  # Bhattacharyya coefficient and Hellinger
  bc <- sum(sqrt(px * py))
  H2 <- max(0, min(1, 1 - bc))   # clip to [0,1]
  if (return_squared) H2 else sqrt(H2)
}
hellinger_discrete <- function(a, b, wa = NULL, wb = NULL,
                               levels_union = NULL,
                               return_squared = FALSE) {
  a <- as.character(a)
  b <- as.character(b)
  if (is.null(levels_union))
    levels_union <- sort(unique(c(a, b)))
  
  # normalize weights within sample
  normalize_w <- function(w, n) {
    if (is.null(w)) rep(1/n, n) else {
      w <- pmax(0, as.numeric(w))
      s <- sum(w)
      if (s == 0) rep(1/n, n) else w / s
    }
  }
  wa <- normalize_w(wa, length(a))
  wb <- normalize_w(wb, length(b))
  
  # empirical probabilities over the same support
  pa <- setNames(numeric(length(levels_union)), levels_union)
  pb <- pa
  ta <- tapply(wa, a, sum)
  tb <- tapply(wb, b, sum)
  pa[names(ta)] <- as.numeric(ta)
  pb[names(tb)] <- as.numeric(tb)
  pa <- pa / sum(pa)
  pb <- pb / sum(pb)
  
  # compute Hellinger
  bc <- sum(sqrt(pa * pb))      # Bhattacharyya coefficient
  H2 <- max(0, min(1, 1 - bc))  # clip to [0,1]
  if (return_squared) H2 else sqrt(H2)
}

# --- Compute Hellinger distances for all covariates and return a table ----
hellinger_summary <- function(nsw_t, X0_matched_original, nsw_c, 
                              covars, cont_vars) {
  discrete_vars <- setdiff(covars, cont_vars)
  results <- data.frame(
    covariate = covars,
    H_matched = NA_real_,
    H_exp = NA_real_,
    type = ifelse(covars %in% cont_vars, "continuous", "discrete"),
    stringsAsFactors = FALSE
  )
  
  for (v in covars) {
    if (v %in% cont_vars) {
      # continuous → histogram Hellinger
      results$H_matched[results$covariate == v] <-
        hellinger_hist(nsw_t[[v]], X0_matched_original[[v]])
      results$H_exp[results$covariate == v] <-
        hellinger_hist(nsw_t[[v]], nsw_c[[v]])
    } else {
      # discrete → discrete Hellinger
      results$H_matched[results$covariate == v] <-
        hellinger_discrete(nsw_t[[v]], X0_matched_original[[v]])
      results$H_exp[results$covariate == v] <-
        hellinger_discrete(nsw_t[[v]], nsw_c[[v]])
    }
  }
  
  # order nicely
  results <- results[order(results$type, results$covariate), ]
  rownames(results) <- NULL
  results
}

sbw_att <- function(dat, ind, out, covars, bal_tol = 0.02, bal_alg = FALSE,
                    bal_std = "group", solver = "quadprog", verbose = FALSE) {
  
  stopifnot(all(c(ind, out, covars) %in% names(dat)))
  
  bal <- list(
    bal_cov = covars,
    bal_alg = bal_alg,       # FALSE = fixed tolerance; TRUE = tuning algorithm
    bal_tol = bal_tol,       # scalar or vector of tolerances
    bal_std = bal_std        # "group" (weighted group SD) or "target"
    # bal_gri, bal_sam keep defaults
  )
  
  fit <- sbw(
    dat = dat,
    ind = ind,
    out = out,
    bal = bal,
    wei = list(wei_sum = TRUE, wei_pos = TRUE),
    sol = list(sol_nam = solver, sol_dis = verbose),
    par = list(par_est = "att", par_tar = NULL),
    mes = verbose
  )
  
  ## Preferred: use package's estimator
  est <- try(estimate(fit), silent = TRUE)
  att_pkg <- if (!inherits(est, "try-error")) {
    if (!is.null(est$att)) est$att else if (!is.null(est$effects$att)) est$effects$att else NA_real_
  } else NA_real_
  
  ## Fallback: compute ATT manually from returned weights
  w <- fit$dat_weights$sbw_weights
  t <- dat[[ind]]; y <- dat[[out]]
  att_manual <- mean(y[t == 1]) - weighted.mean(y[t == 0], w[t == 0])
  
  list(att = if (is.finite(att_pkg)) att_pkg else att_manual,
       att_pkg = att_pkg,
       att_manual = att_manual,
       fit = fit)
}
