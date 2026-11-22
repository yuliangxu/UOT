#!/usr/bin/env Rscript

# Figure 1, from the repository root:
#   Rscript --vanilla simulation/ET_round2/MMD_heatmap.R
#   Rscript --vanilla simulation/ET_round2/MMD_heatmap.R --smoke-test
# --reps=N runs or resumes the first N original seeds in each grid cell.

options(stringsAsFactors = FALSE)

parse_arguments <- function(args) {
  value_after_equals <- function(prefix, default) {
    hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
    if (!length(hit)) return(default)
    sub(paste0("^", prefix, "="), "", hit[[length(hit)]])
  }

  list(
    smoke_test = "--smoke-test" %in% args,
    validate_mmd = !"--skip-validation" %in% args,
    reps = as.integer(value_after_equals("--reps", "50")),
    seed = as.integer(value_after_equals("--seed", "20260813")),
    output_dir = value_after_equals(
      "--output-dir", "simulation/ET_round2/figure1"
    ),
    plot_dir = value_after_equals(
      "--plot-dir", "simulation/ET_round2/figure1"
    )
  )
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing R package(s): ", paste(missing, collapse = ", "),
      ". Install them before running this experiment.",
      call. = FALSE
    )
  }
}

assert_repository_root <- function() {
  expected <- file.path("simulation", "ET_round2", "MMD_heatmap.R")
  if (!file.exists(expected)) {
    stop("Run this script from the repository root. Expected to find ", expected,
         call. = FALSE)
  }
}

# Marginals have bounded eigenvalues, different covariance structures, and a
# dimension-independent Mahalanobis-scale mean separation. This gives them
# substantial overlapping support for every d in the experiment.
make_overlapping_marginals <- function(d) {
  stopifnot(length(d) == 1L, d >= 2L)

  index <- seq_len(d)
  A <- outer(index, index, function(i, j) 0.30^abs(i - j))
  B <- outer(index, index, function(i, j) 0.18^abs(i - j))

  scale_a <- seq(0.90, 1.10, length.out = d)
  scale_b <- rev(seq(0.85, 1.15, length.out = d))
  A <- diag(sqrt(scale_a), d) %*% A %*% diag(sqrt(scale_a), d)
  B <- diag(sqrt(scale_b), d) %*% B %*% diag(sqrt(scale_b), d)

  # A small deterministic rotation makes B distinct without degrading its
  # conditioning. Apply rotations to adjacent coordinate pairs.
  Q <- diag(d)
  theta <- pi / 12
  for (i in seq(1L, d - 1L, by = 2L)) {
    Q[c(i, i + 1L), c(i, i + 1L)] <- matrix(
      c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2L, 2L
    )
  }
  B <- Q %*% B %*% t(Q)

  direction <- rep(c(1, -1), length.out = d)
  direction <- direction / sqrt(sum(direction^2))
  mean_difference <- 0.60 * direction

  list(
    a = -mean_difference / 2,
    b = mean_difference / 2,
    A = (A + t(A)) / 2,
    B = (B + t(B)) / 2
  )
}

# Use the same Janati et al. Gaussian benchmark as the Appendix C experiment.
gaussian <- new.env(parent = environment())
sys.source("simulation/ET_round2/CLT_gaussian_experiment.R", envir = gaussian)
squared_distances <- gaussian$squared_distances

# G has X1 in its rows and X0 in its columns, as required by the MMD routine.
compute_uot_plan <- function(X0, X1, epsilon, rho, max_iter = 2000L) {
  cost <- squared_distances(X1, X0)
  log_a <- -log(nrow(X1))
  log_b <- -log(nrow(X0))
  tau <- rho / (rho + epsilon)
  log_kernel <- -cost / epsilon
  row_lse <- function(x) {
    peak <- apply(x, 1L, max)
    peak + log(rowSums(exp(x - peak)))
  }

  log_v <- numeric(nrow(X0))
  for (iteration in seq_len(max_iter)) {
    log_u <- -tau * (row_lse(sweep(log_kernel, 2L, log_v, "+")) + log_b)
    updated_v <- -tau * (row_lse(t(log_kernel + log_u)) + log_a)
    step <- max(abs(updated_v - log_v))
    log_v <- updated_v
    if (epsilon * step > 1e-12) next

    log_plan <- sweep(log_kernel + log_u, 2L, log_v, "+") + log_a + log_b
    log_rows <- row_lse(log_plan)
    log_cols <- row_lse(t(log_plan))
    # First-order derivative for every raw mass, including cells below the
    # floating-point range. Normalize only after checking the original problem.
    derivative <- cost + epsilon * (log_plan - log_a - log_b) +
      rho * outer(log_rows - log_a, log_cols - log_b, "+")
    residual <- max(abs(derivative))
    if (!is.finite(residual)) stop("Nonfinite KL first-order residual.")
    if (residual > 1e-10) next
    peak <- max(log_rows)
    log_mass <- peak + log(sum(exp(log_rows - peak)))
    G <- exp(log_plan - log_mass)
    if (!is.finite(log_mass) || log_mass > 1e-10 || abs(sum(G) - 1) > 1e-10)
      stop("Invalid raw mass or normalized UOT plan.")
    return(list(G = G, max_abs_foc = residual, log_mass = log_mass, iterations = iteration))
  }
  stop("Log-domain Sinkhorn did not meet the first-order tolerance.")
}

# Exact squared MMD. For z_ij=(x0_i,x1_j), the Gaussian kernel factorizes:
# k(z_ij,z_kl)=K0[i,k] K1[j,l]. Thus the empirical-empirical term is
# sum(G * (K1 %*% G %*% K0)) rather than an (n^2)-by-(n^2) calculation.
factorized_closed_form_mmd <- function(X0, X1, G, mu, Sigma) {
  d <- ncol(X0)
  joint_d <- 2L * d

  K0 <- exp(-squared_distances(X0) / joint_d)
  K1 <- exp(-squared_distances(X1) / joint_d)
  E_XX <- sum(G * (K1 %*% G %*% K0))

  gaussian_kernel_cov <- diag(joint_d) + 2 * Sigma / joint_d
  gaussian_kernel_cov <- (gaussian_kernel_cov + t(gaussian_kernel_cov)) / 2
  chol_cov <- chol(gaussian_kernel_cov)
  cov_inv <- chol2inv(chol_cov)
  log_det_cov <- 2 * sum(log(diag(chol_cov)))

  # Evaluate the Gaussian cross-term at every transport support pair without
  # constructing any pairwise (n^2)-by-(n^2) matrix.
  pairs <- cbind(
    X0[rep(seq_len(nrow(X0)), each = nrow(X1)), , drop = FALSE],
    X1[rep(seq_len(nrow(X1)), times = nrow(X0)), , drop = FALSE]
  )
  centered <- sweep(pairs, 2L, mu, "-")
  quadratic <- rowSums((centered %*% cov_inv) * centered)
  pair_weights <- as.vector(G)
  E_XY <- sum(pair_weights * exp(-quadratic / joint_d - log_det_cov / 2))

  # If Y,Y' are iid N(mu,Sigma), Y-Y' is N(0,2 Sigma).
  yy_cov <- diag(joint_d) + 4 * Sigma / joint_d
  yy_chol <- chol((yy_cov + t(yy_cov)) / 2)
  E_YY <- exp(-sum(log(diag(yy_chol))))

  mmd2 <- E_XX - 2 * E_XY + E_YY
  list(mmd2 = mmd2, E_XX = E_XX, E_XY = E_XY, E_YY = E_YY)
}

expanded_empirical_term <- function(X0, X1, G) {
  d <- ncol(X0)
  pairs <- cbind(
    X0[rep(seq_len(nrow(X0)), each = nrow(X1)), , drop = FALSE],
    X1[rep(seq_len(nrow(X1)), times = nrow(X0)), , drop = FALSE]
  )
  weights <- as.vector(G)
  kernel <- exp(-squared_distances(pairs) / (2L * d))
  sum(outer(weights, weights) * kernel)
}

validate_factorization <- function() {
  set.seed(17)
  X0 <- matrix(rnorm(15), 5, 3)
  X1 <- matrix(rnorm(12), 4, 3)
  G <- matrix(rexp(20), 4, 5)
  G <- G / sum(G)
  K0 <- exp(-squared_distances(X0) / 6)
  K1 <- exp(-squared_distances(X1) / 6)
  factorized <- sum(G * (K1 %*% G %*% K0))
  expanded <- expanded_empirical_term(X0, X1, G)
  if (!isTRUE(all.equal(factorized, expanded, tolerance = 1e-12))) {
    stop("Factorized MMD validation failed.", call. = FALSE)
  }
  message("Factorized MMD validation passed (absolute error: ",
          format(abs(factorized - expanded), scientific = TRUE), ").")
}

cell_seed <- function(master_seed, d, n, replication) {
  # Stable across grid order and resumptions; kept below R's integer limit.
  as.integer((master_seed + 1000003 * d + 1009 * n + replication) %% 2147483646 + 1)
}

run_one_replication <- function(d, n, replication, master_seed, epsilon, rho) {
  seed <- cell_seed(master_seed, d, n, replication)
  set.seed(seed)
  marginal <- make_overlapping_marginals(d)
  truth <- gaussian$gaussian_uot_population(
    marginal$a, marginal$b, marginal$A, marginal$B, epsilon, rho
  )
  X0 <- MASS::mvrnorm(n, mu = marginal$a, Sigma = marginal$A)
  X1 <- MASS::mvrnorm(n, mu = marginal$b, Sigma = marginal$B)

  start_total <- proc.time()[["elapsed"]]
  start_sinkhorn <- proc.time()[["elapsed"]]
  fit <- compute_uot_plan(X0, X1, epsilon, rho)
  sinkhorn_seconds <- proc.time()[["elapsed"]] - start_sinkhorn

  start_mmd <- proc.time()[["elapsed"]]
  mmd <- factorized_closed_form_mmd(X0, X1, fit$G, truth$mean, truth$covariance)
  mmd_seconds <- proc.time()[["elapsed"]] - start_mmd

  if (!is.finite(mmd$mmd2) || mmd$mmd2 < -1e-12) stop("Invalid squared MMD.")
  data.frame(
    d = d,
    n = n,
    replication = replication,
    seed = seed,
    max_abs_foc = fit$max_abs_foc,
    log_mass = fit$log_mass,
    iterations = fit$iterations,
    largest_pair_weight = max(fit$G),
    mmd2 = mmd$mmd2,
    mmd = sqrt(max(mmd$mmd2, 0)),
    E_XX = mmd$E_XX,
    E_XY = mmd$E_XY,
    E_YY = mmd$E_YY,
    sinkhorn_seconds = sinkhorn_seconds,
    mmd_seconds = mmd_seconds,
    total_seconds = proc.time()[["elapsed"]] - start_total
  )
}

atomic_save_rds <- function(object, path) {
  temporary <- tempfile(pattern = "checkpoint-", tmpdir = dirname(path), fileext = ".rds")
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not replace checkpoint: ", path, call. = FALSE)
  }
}

summarize_results <- function(results) {
  keys <- interaction(results$d, results$n, drop = TRUE)
  pieces <- lapply(split(results, keys), function(x) {
    data.frame(
      d = x$d[[1]],
      n = x$n[[1]],
      replications = nrow(x),
      median_mmd2 = median(x$mmd2),
      q25_mmd2 = unname(quantile(x$mmd2, 0.25)),
      q75_mmd2 = unname(quantile(x$mmd2, 0.75)),
      median_mmd = median(x$mmd),
      q25_mmd = unname(quantile(x$mmd, 0.25)),
      q75_mmd = unname(quantile(x$mmd, 0.75)),
      median_seconds = median(x$total_seconds)
    )
  })
  summary <- do.call(rbind, pieces)
  summary[order(summary$d, summary$n), , drop = FALSE]
}

make_heatmap <- function(summary, path, target_reps) {
  summary$d_factor <- factor(summary$d, levels = sort(unique(summary$d)))
  summary$n_factor <- factor(summary$n, levels = sort(unique(summary$n)))

  plot <- ggplot2::ggplot(
    summary, ggplot2::aes(x = n_factor, y = d_factor, fill = median_mmd)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", median_mmd),
                   colour = median_mmd > mean(range(median_mmd))), size = 4.2
    ) +
    ggplot2::scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
    ggplot2::scale_fill_viridis_c(option = "C", direction = -1) +
    ggplot2::labs(
      title = paste0("Median MMD over ", target_reps, " Replications"),
      subtitle = expression(epsilon == 0.01 * "," ~~ rho == 0.002),
      x = "Sample size per marginal (n)",
      y = "Dimension (d)",
      fill = "Median MMD"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 18, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 15),
      axis.title = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 13)
    )

  ggplot2::ggsave(path, plot, width = 8.5, height = 5.5, units = "in")
  plot
}

make_lineplot <- function(summary, path, target_reps) {
  summary$d_factor <- factor(summary$d, levels = sort(unique(summary$d)))

  plot <- ggplot2::ggplot(
    summary,
    ggplot2::aes(
      x = sqrt(n), y = median_mmd, color = d_factor, group = d_factor
    )
  ) +
    ggplot2::geom_line(linewidth = 1.0) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = q25_mmd, ymax = q75_mmd),
      width = 0.12, linewidth = 0.8
    ) +
    ggplot2::scale_x_continuous(
      breaks = sqrt(sort(unique(summary$n))),
      labels = function(x) format(round(x, 2), trim = TRUE)
    ) +
    ggplot2::labs(
      title = paste0(
        "Median MMD with 25% and 75% Quantiles over ",
        target_reps, " Replications"
      ),
      subtitle = expression(epsilon == 0.01 * "," ~~ rho == 0.002),
      x = expression("Square-root sample size per marginal " * sqrt(n)),
      y = "MMD",
      color = "d"
    ) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(size = 18, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 15),
      axis.title = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 13),
      legend.title = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 13)
    )

  ggplot2::ggsave(path, plot, width = 8.5, height = 5.5, units = "in")
  plot
}

main <- function() {
  args <- parse_arguments(commandArgs(trailingOnly = TRUE))
  assert_repository_root()
  require_packages(c("MASS", "ggplot2"))

  if (is.na(args$reps) || args$reps < 1L) {
    stop("--reps must be a positive integer.", call. = FALSE)
  }
  if (args$validate_mmd) validate_factorization()

  epsilon <- 0.01
  rho <- 0.002
  if (args$smoke_test) {
    d_grid <- c(2L, 4L)
    n_grid <- c(10L, 15L)
    target_reps <- 1L
    run_name <- "smoke_test"
  } else {
    d_grid <- c(2L, 3L, 4L, 5L, 6L, 8L, 10L)
    n_grid <- c(20L, 40L, 60L, 80L, 100L, 150L, 200L)
    target_reps <- args$reps
    run_name <- "mmd_heatmap"
  }

  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(args$plot_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint_path <- file.path(args$output_dir, paste0(run_name, "_checkpoint.rds"))
  summary_path <- file.path(args$output_dir, paste0(run_name, "_summary.csv"))
  metadata_path <- file.path(args$output_dir, paste0(run_name, "_metadata.rds"))
  plot_path <- file.path(args$plot_dir, paste0(run_name, ".pdf"))
  lineplot_path <- file.path(args$plot_dir, paste0(run_name, "_lineplot.pdf"))

  if (file.exists(checkpoint_path) && !file.exists(metadata_path))
    stop("Checkpoint has no matching specification. Choose a fresh output directory.")
  if (file.exists(checkpoint_path)) {
    results <- readRDS(checkpoint_path)
  } else {
    results <- data.frame()
  }

  metadata <- list(
    d_grid = d_grid,
    n_grid = n_grid,
    target_replications = target_reps,
    specification = "squared-euclidean-product-kl-janati-v1",
    master_seed = args$seed,
    epsilon = epsilon,
    rho = rho,
    kernel = "exp(-squared_distance/(2*d)) on the 2d-dimensional coupling",
    marginal_specification = "distinct well-conditioned overlapping Gaussians",
    created_at = Sys.time(),
    R_version = R.version.string
  )
  if (file.exists(metadata_path)) {
    previous_metadata <- readRDS(metadata_path)
    compatible <- identical(previous_metadata$specification, metadata$specification) &&
      identical(previous_metadata$d_grid, metadata$d_grid) &&
      identical(previous_metadata$n_grid, metadata$n_grid) &&
      identical(previous_metadata$master_seed, metadata$master_seed) &&
      identical(previous_metadata$epsilon, metadata$epsilon) &&
      identical(previous_metadata$rho, metadata$rho)
    if (!compatible) {
      stop(
        "The existing checkpoint uses a different cost, solver, grid, seed, ",
        "epsilon, or rho. Choose a different --output-dir or remove the old ",
        "smoke-test output before changing its configuration.", call. = FALSE
      )
    }
    metadata$created_at <- previous_metadata$created_at
  }
  atomic_save_rds(metadata, metadata_path)

  design <- expand.grid(
    d = d_grid, n = n_grid, replication = seq_len(target_reps),
    KEEP.OUT.ATTRS = FALSE
  )
  if (nrow(results)) {
    if (!all(c("seed", "max_abs_foc", "largest_pair_weight") %in% names(results)))
      stop("Checkpoint lacks the Figure 1 verification columns. Use a fresh output directory.")
    completed <- paste(results$d, results$n, results$replication, sep = ":")
    design_key <- paste(design$d, design$n, design$replication, sep = ":")
    design <- design[!design_key %in% completed, , drop = FALSE]
  }

  message("Completed: ", nrow(results), "; remaining: ", nrow(design), ".")
  # Save after each replication ID (at most 49 jobs) to resume long runs.
  for (jobs in split(design, design$replication)) {
    new_rows <- lapply(seq_len(nrow(jobs)), function(i) {
      run_one_replication(jobs$d[i], jobs$n[i], jobs$replication[i],
                          args$seed, epsilon, rho)
    })
    results <- rbind(results, do.call(rbind, new_rows))
    atomic_save_rds(results, checkpoint_path)
    message("Completed ", nrow(results), " fits.")
  }
  selected <- results[results$d %in% d_grid & results$n %in% n_grid &
                      results$replication <= target_reps, , drop = FALSE]
  selected <- selected[order(selected$d, selected$n, selected$replication), , drop = FALSE]
  stopifnot(nrow(selected) == length(d_grid) * length(n_grid) * target_reps,
            !anyDuplicated(selected[c("d", "n", "replication")]),
            all(is.finite(selected$mmd)), all(selected$max_abs_foc <= 1e-10))
  summary <- summarize_results(selected)
  utils::write.csv(selected, file.path(args$output_dir, paste0(run_name, "_replications.csv")),
                   row.names = FALSE)
  utils::write.csv(summary, summary_path, row.names = FALSE)
  make_heatmap(summary, plot_path, target_reps)
  make_lineplot(summary, lineplot_path, target_reps)
  message("Saved checkpoint: ", checkpoint_path)
  message("Saved summary: ", summary_path)
  message("Saved heatmap: ", plot_path)
  message("Saved line plot: ", lineplot_path)
}

if (sys.nframe() == 0L) main()
