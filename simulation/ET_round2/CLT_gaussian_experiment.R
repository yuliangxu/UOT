#!/usr/bin/env Rscript

# ET round 2: numerical illustration of the Gaussian-process limit in Theorem 1.
#
# Run from the repository root. The Python selected by reticulate must contain
# POT (the `ot` module):
#
#   export RETICULATE_PYTHON=/usr/bin/python3
#   Rscript simulation/ET_round2/CLT_gaussian_experiment.R --smoke-test
#   Rscript simulation/ET_round2/CLT_gaussian_experiment.R --reps=100
#   Rscript simulation/ET_round2/CLT_gaussian_experiment.R --reps=2000
#
# Increasing --reps resumes both the evaluation checkpoint and a separately
# seeded reference checkpoint at the largest n. Thus --reps=2000 gives 2,000
# evaluation runs at every n plus 2,000 independent large-n reference runs.
# The experiment uses the unnormalized UOT plan, fixed squared Euclidean cost,
# epsilon=0.5, rho=1, and the KL reference alpha (x) beta, matching
# arXiv:2006.02572.

options(stringsAsFactors = FALSE)

parse_arguments <- function(args) {
  value_after_equals <- function(prefix, default) {
    hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
    if (!length(hit)) return(default)
    sub(paste0("^", prefix, "="), "", hit[[length(hit)]])
  }

  list(
    smoke_test = "--smoke-test" %in% args,
    reps = as.integer(value_after_equals("--reps", "2000")),
    seed = as.integer(value_after_equals("--seed", "20260814")),
    output_dir = value_after_equals(
      "--output-dir", "simulation/ET_round2/result"
    ),
    plot_dir = value_after_equals(
      "--plot-dir", "simulation/ET_round2/plot"
    )
  )
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing R package(s): ", paste(missing, collapse = ", "),
      ". Install them before running this experiment.", call. = FALSE
    )
  }
}

assert_repository_root <- function() {
  expected <- file.path("simulation", "ET_round2", "CLT_gaussian_experiment.R")
  if (!file.exists(expected)) {
    stop("Run this script from the repository root. Expected to find ", expected,
         call. = FALSE)
  }
}

squared_distances <- function(X, Y = X) {
  D2 <- outer(rowSums(X^2), rowSums(Y^2), "+") - 2 * tcrossprod(X, Y)
  pmax(D2, 0)
}

principal_matrix_sqrt <- function(M, tolerance = 1e-8) {
  decomposition <- eigen(M)
  if (max(abs(Im(decomposition$values))) > tolerance ||
      min(Re(decomposition$values)) < -tolerance) {
    stop("The requested principal matrix square root is not real and positive.",
         call. = FALSE)
  }
  root <- decomposition$vectors %*%
    diag(sqrt(pmax(Re(decomposition$values), 0)), nrow(M)) %*%
    solve(decomposition$vectors)
  root <- Re(root)
  relative_residual <- norm(root %*% root - M, type = "F") /
    max(1, norm(M, type = "F"))
  if (!is.finite(relative_residual) || relative_residual > 1e-7) {
    stop("The principal matrix square root failed its residual check.",
         call. = FALSE)
  }
  root
}

log_positive_determinant <- function(M, label) {
  determinant_result <- determinant(M, logarithm = TRUE)
  if (determinant_result$sign <= 0) {
    stop(label, " does not have a positive determinant.", call. = FALSE)
  }
  as.numeric(determinant_result$modulus)
}

# Closed form from Janati et al. (2020), arXiv:2006.02572. Their entropy
# coefficient is 2*sigma^2, which is called epsilon here; gamma is called rho.
gaussian_uot_population <- function(a, b, A, B, epsilon, rho,
                                    mass_a = 1, mass_b = 1) {
  d <- length(a)
  stopifnot(
    length(b) == d, all(dim(A) == c(d, d)), all(dim(B) == c(d, d)),
    epsilon > 0, rho > 0, mass_a > 0, mass_b > 0
  )

  sigma2 <- epsilon / 2
  sigma <- sqrt(sigma2)
  gamma <- rho
  tau <- gamma / (2 * sigma2 + gamma)
  lambda <- sigma2 + gamma / 2
  identity <- diag(d)

  X <- A + B + lambda * identity
  X_inv <- chol2inv(chol((X + t(X)) / 2))
  A_tilde <- gamma / 2 * (
    identity - lambda * chol2inv(chol(A + lambda * identity))
  )
  B_tilde <- gamma / 2 * (
    identity - lambda * chol2inv(chol(B + lambda * identity))
  )
  C <- principal_matrix_sqrt(
    A_tilde %*% B_tilde / tau + sigma2^2 / 4 * identity
  ) - sigma2 / 2 * identity

  mean <- c(
    a + A %*% X_inv %*% (b - a),
    b - B %*% X_inv %*% (b - a)
  )

  identity_C <- identity + C / lambda
  identity_C_transpose <- identity + t(C) / lambda
  covariance <- rbind(
    cbind(
      identity_C %*% (A - A %*% X_inv %*% A),
      C + identity_C %*% A %*% X_inv %*% B
    ),
    cbind(
      t(C) + identity_C_transpose %*% B %*% X_inv %*% A,
      identity_C_transpose %*% (B - B %*% X_inv %*% B)
    )
  )
  asymmetry <- norm(covariance - t(covariance), type = "F") /
    max(1, norm(covariance, type = "F"))
  if (!is.finite(asymmetry) || asymmetry > 1e-7) {
    stop("The Gaussian UOT covariance is not numerically symmetric.", call. = FALSE)
  }
  covariance <- (Re(covariance) + t(Re(covariance))) / 2
  if (min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("The Gaussian UOT covariance is not positive definite.", call. = FALSE)
  }

  mean_difference <- b - a
  mean_quadratic <- as.numeric(t(mean_difference) %*% X_inv %*% mean_difference)
  log_mass <- d * sigma2 / (gamma + sigma2) * log(sigma) +
    (
      log(mass_a * mass_b) +
        log_positive_determinant(C, "C") +
        0.5 * (
          tau * log_positive_determinant(A_tilde %*% B_tilde, "A_tilde B_tilde") -
            log_positive_determinant(A %*% B, "A B")
        )
    ) / (tau + 1) -
    mean_quadratic / (2 * (tau + 1)) -
    0.5 * log_positive_determinant(
      C - 2 / gamma * A_tilde %*% B_tilde,
      "C - 2 A_tilde B_tilde / gamma"
    )

  list(
    mass = exp(log_mass),
    log_mass = log_mass,
    mean = as.numeric(mean),
    covariance = covariance,
    C = C,
    tau = tau,
    lambda = lambda,
    sigma2 = sigma2
  )
}

make_marginals <- function() {
  list(
    a = c(-0.4, 0.1),
    b = c(0.5, -0.3),
    A = matrix(c(0.7, -0.5, -0.5, 1.0), nrow = 2),
    B = matrix(c(1.0, -0.34, -0.34, 0.7), nrow = 2)
  )
}

make_test_functions <- function(population, grid_size = 7L) {
  d <- length(population$mean) / 2L
  projection_spec <- data.frame(
    projection = c("p11", "p12", "p21", "p22"),
    source_coordinate = c(1L, 1L, 2L, 2L),
    target_coordinate = c(1L, 2L, 1L, 2L),
    projection_label = c(
      "(1,1): matched", "(1,2): cross",
      "(2,1): cross", "(2,2): matched"
    )
  )

  projections <- lapply(seq_len(nrow(projection_spec)), function(k) {
    specification <- projection_spec[k, ]
    projection_indices <- c(
      specification$source_coordinate,
      d + specification$target_coordinate
    )
    projected_mean <- population$mean[projection_indices]
    projected_covariance <- population$covariance[
      projection_indices, projection_indices, drop = FALSE
    ]
    projected_sd <- sqrt(diag(projected_covariance))
    source_centers <- seq(
      projected_mean[[1]] - 2 * projected_sd[[1]],
      projected_mean[[1]] + 2 * projected_sd[[1]],
      length.out = grid_size
    )
    target_centers <- seq(
      projected_mean[[2]] - 2 * projected_sd[[2]],
      projected_mean[[2]] + 2 * projected_sd[[2]],
      length.out = grid_size
    )
    bandwidth <- 0.5 * mean(projected_sd)
    grid <- expand.grid(
      source_center = source_centers,
      target_center = target_centers,
      KEEP.OUT.ATTRS = FALSE
    )
    grid$source_index <- match(grid$source_center, source_centers)
    grid$target_index <- match(grid$target_center, target_centers)
    grid$projection <- specification$projection
    grid$projection_label <- specification$projection_label
    grid$source_coordinate <- specification$source_coordinate
    grid$target_coordinate <- specification$target_coordinate
    grid$feature <- sprintf(
      "%s_bump_s%02d_t%02d", specification$projection,
      grid$source_index, grid$target_index
    )

    kernel_covariance <- projected_covariance + bandwidth^2 * diag(2)
    kernel_covariance_inv <- chol2inv(chol(kernel_covariance))
    determinant_factor <- exp(
      -0.5 * log_positive_determinant(
        diag(2) + projected_covariance / bandwidth^2,
        paste0("I + projected covariance / bandwidth^2 for ",
               specification$projection)
      )
    )
    population_bumps <- apply(
      grid[, c("source_center", "target_center"), drop = FALSE],
      1,
      function(center) {
        difference <- center - projected_mean
        population$mass * determinant_factor *
          exp(-0.5 * as.numeric(
            t(difference) %*% kernel_covariance_inv %*% difference
          ))
      }
    )

    list(
      specification = specification,
      projection_indices = projection_indices,
      projected_mean = projected_mean,
      projected_covariance = projected_covariance,
      source_centers = source_centers,
      target_centers = target_centers,
      bandwidth = bandwidth,
      grid = grid,
      population_integrals = population_bumps
    )
  })
  names(projections) <- projection_spec$projection
  grid <- do.call(rbind, lapply(projections, `[[`, "grid"))
  rownames(grid) <- NULL
  population_bumps <- unlist(
    lapply(projections, `[[`, "population_integrals"), use.names = FALSE
  )
  feature_names <- c("total_mass", grid$feature)
  population_integrals <- c(total_mass = population$mass, population_bumps)
  names(population_integrals) <- feature_names

  list(
    grid = grid,
    projection_spec = projection_spec,
    projections = projections,
    population_integrals = population_integrals,
    feature_names = feature_names
  )
}

compute_empirical_integrals <- function(X0, X1, plan, test_functions) {
  bump_integrals <- unlist(lapply(test_functions$projections, function(projection) {
    source_coordinate <- projection$specification$source_coordinate
    target_coordinate <- projection$specification$target_coordinate
    source_values <- outer(
      X0[, source_coordinate], projection$source_centers,
      function(x, center) {
        exp(-(x - center)^2 / (2 * projection$bandwidth^2))
      }
    )
    target_values <- outer(
      X1[, target_coordinate], projection$target_centers,
      function(x, center) {
        exp(-(x - center)^2 / (2 * projection$bandwidth^2))
      }
    )
    as.vector(t(source_values) %*% plan %*% target_values)
  }), use.names = FALSE)
  values <- c(total_mass = sum(plan), bump_integrals)
  names(values) <- test_functions$feature_names
  values
}

compute_uot_plan <- function(X0, X1, epsilon, rho, max_iter = 10000L,
                             tolerance = 1e-10) {
  n0 <- nrow(X0)
  n1 <- nrow(X1)
  cost <- squared_distances(X0, X1)
  weights0 <- rep(1 / n0, n0)
  weights1 <- rep(1 / n1, n1)
  ot <- reticulate::import("ot", convert = TRUE)

  result <- ot$unbalanced$sinkhorn_unbalanced(
    weights0, weights1, cost,
    reg = epsilon,
    reg_m = rho,
    method = "sinkhorn_stabilized",
    reg_type = "kl",
    numItermax = as.integer(max_iter),
    stopThr = tolerance,
    log = TRUE
  )
  plan <- as.matrix(result[[1]])
  log <- result[[2]]
  errors <- unlist(log$err)
  final_error <- if (length(errors)) tail(errors, 1) else NA_real_
  converged <- is.finite(final_error) && final_error <= tolerance * 10

  if (!all(dim(plan) == c(n0, n1)) || !all(is.finite(plan)) ||
      any(plan < 0) || sum(plan) <= 0) {
    stop("POT returned an invalid unbalanced transport plan.", call. = FALSE)
  }

  list(
    plan = plan,
    iterations = length(errors),
    final_error = final_error,
    converged = converged
  )
}

cell_seed <- function(master_seed, n, replication, seed_stream = 0L) {
  # Stream 0 reproduces the original evaluation seeds. Stream 1 is reserved
  # for the independent large-n covariance-reference batch.
  as.integer(
    (master_seed + 1000003 * n + replication + 1000000007 * seed_stream) %%
      2147483646 + 1
  )
}

run_one_replication <- function(n, replication, master_seed, marginals,
                                population, test_functions, epsilon, rho,
                                seed_stream = 0L) {
  seed <- cell_seed(master_seed, n, replication, seed_stream)
  set.seed(seed)
  X0 <- MASS::mvrnorm(n, mu = marginals$a, Sigma = marginals$A)
  X1 <- MASS::mvrnorm(n, mu = marginals$b, Sigma = marginals$B)

  start <- proc.time()[["elapsed"]]
  solution <- compute_uot_plan(X0, X1, epsilon, rho)
  integrals <- compute_empirical_integrals(
    X0, X1, solution$plan, test_functions
  )
  elapsed <- proc.time()[["elapsed"]] - start

  diagnostics <- data.frame(
    n = n,
    replication = replication,
    seed = seed,
    plan_mass = sum(solution$plan),
    iterations = solution$iterations,
    final_error = solution$final_error,
    converged = solution$converged,
    elapsed_seconds = elapsed
  )
  list(diagnostics = diagnostics, integrals = integrals)
}

atomic_save_rds <- function(object, path) {
  temporary <- tempfile(pattern = "checkpoint-", tmpdir = dirname(path), fileext = ".rds")
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not replace checkpoint: ", path, call. = FALSE)
  }
}

empty_checkpoint <- function(feature_names) {
  list(
    diagnostics = data.frame(
      n = integer(), replication = integer(), seed = integer(),
      plan_mass = double(), iterations = integer(), final_error = double(),
      converged = logical(), elapsed_seconds = double()
    ),
    integrals = matrix(
      numeric(), nrow = 0, ncol = length(feature_names),
      dimnames = list(NULL, feature_names)
    )
  )
}

validate_checkpoint <- function(checkpoint, feature_names, allowed_n,
                                master_seed, seed_stream, label) {
  if (!is.list(checkpoint) || is.null(checkpoint$diagnostics) ||
      is.null(checkpoint$integrals)) {
    stop(label, " checkpoint has an invalid structure.", call. = FALSE)
  }
  if (nrow(checkpoint$diagnostics) != nrow(checkpoint$integrals)) {
    stop(label, " checkpoint diagnostics and integrals are misaligned.",
         call. = FALSE)
  }
  if (!identical(colnames(checkpoint$integrals), feature_names)) {
    stop(label, " checkpoint test-function columns do not match.",
         call. = FALSE)
  }
  if (!nrow(checkpoint$diagnostics)) return(invisible(TRUE))
  if (!all(checkpoint$diagnostics$n %in% allowed_n)) {
    stop(label, " checkpoint contains an unexpected sample size.",
         call. = FALSE)
  }
  keys <- paste(
    checkpoint$diagnostics$n, checkpoint$diagnostics$replication, sep = ":"
  )
  if (anyDuplicated(keys) || anyDuplicated(checkpoint$diagnostics$seed)) {
    stop(label, " checkpoint contains duplicate jobs or seeds.", call. = FALSE)
  }
  expected_seeds <- mapply(
    cell_seed,
    n = checkpoint$diagnostics$n,
    replication = checkpoint$diagnostics$replication,
    MoreArgs = list(master_seed = master_seed, seed_stream = seed_stream)
  )
  if (!identical(as.integer(checkpoint$diagnostics$seed),
                 as.integer(expected_seeds))) {
    stop(label, " checkpoint does not match its configured seed stream.",
         call. = FALSE)
  }
  invisible(TRUE)
}

append_result <- function(checkpoint, result) {
  checkpoint$diagnostics <- rbind(checkpoint$diagnostics, result$diagnostics)
  checkpoint$integrals <- rbind(checkpoint$integrals, result$integrals)
  ordering <- order(
    checkpoint$diagnostics$n, checkpoint$diagnostics$replication
  )
  checkpoint$diagnostics <- checkpoint$diagnostics[ordering, , drop = FALSE]
  checkpoint$integrals <- checkpoint$integrals[ordering, , drop = FALSE]
  checkpoint
}

complete_design <- function(checkpoint, design, checkpoint_path, batch_label,
                            seed_stream, master_seed, marginals, population,
                            test_functions, epsilon, rho) {
  message(
    batch_label, " completed: ", nrow(checkpoint$diagnostics),
    "; remaining: ", nrow(design), "."
  )
  for (row in seq_len(nrow(design))) {
    job <- design[row, ]
    if (row == 1L || row %% 25L == 0L || row == nrow(design)) {
      message(
        "[", batch_label, " ", row, "/", nrow(design), "] n=", job$n,
        ", replication=", job$replication
      )
    }
    result <- run_one_replication(
      job$n, job$replication, master_seed, marginals, population,
      test_functions, epsilon, rho, seed_stream
    )
    checkpoint <- append_result(checkpoint, result)
    if (row %% 10L == 0L || row == nrow(design)) {
      atomic_save_rds(checkpoint, checkpoint_path)
    }
  }
  checkpoint
}

compute_fluctuations <- function(checkpoint, population_integrals) {
  centered <- sweep(checkpoint$integrals, 2, population_integrals, "-")
  centered * sqrt(checkpoint$diagnostics$n)
}

summarize_replications <- function(checkpoint, fluctuations, population_mass) {
  pieces <- lapply(split(seq_len(nrow(fluctuations)), checkpoint$diagnostics$n), function(i) {
    data.frame(
      n = checkpoint$diagnostics$n[i[[1]]],
      replications = length(i),
      mean_plan_mass = mean(checkpoint$diagnostics$plan_mass[i]),
      population_mass = population_mass,
      mass_bias = mean(checkpoint$diagnostics$plan_mass[i]) - population_mass,
      median_seconds = median(checkpoint$diagnostics$elapsed_seconds[i]),
      convergence_rate = mean(checkpoint$diagnostics$converged[i]),
      max_absolute_mean_fluctuation = max(abs(colMeans(fluctuations[i, , drop = FALSE])))
    )
  })
  do.call(rbind, pieces)
}

heatmap_data <- function(values, test_functions, n_values, statistic) {
  pieces <- lapply(seq_along(n_values), function(k) {
    data.frame(
      test_functions$grid,
      n = n_values[[k]],
      value = values[k, ],
      statistic = statistic
    )
  })
  do.call(rbind, pieces)
}

make_fluctuation_heatmaps <- function(checkpoint, fluctuations, test_functions,
                                      plot_dir, prefix) {
  n_values <- sort(unique(checkpoint$diagnostics$n))
  bump_columns <- seq_len(nrow(test_functions$grid)) + 1L
  mean_values <- t(vapply(n_values, function(n) {
    indices <- which(checkpoint$diagnostics$n == n)
    colMeans(fluctuations[indices, bump_columns, drop = FALSE])
  }, numeric(length(bump_columns))))
  sd_values <- t(vapply(n_values, function(n) {
    indices <- which(checkpoint$diagnostics$n == n)
    apply(fluctuations[indices, bump_columns, drop = FALSE], 2, stats::sd)
  }, numeric(length(bump_columns))))
  representative_values <- t(vapply(n_values, function(n) {
    index <- which(checkpoint$diagnostics$n == n)[[1]]
    fluctuations[index, bump_columns]
  }, numeric(length(bump_columns))))

  plot_heatmap <- function(data, title, fill_label, diverging, filename) {
    data$projection_label <- factor(
      data$projection_label,
      levels = test_functions$projection_spec$projection_label
    )
    plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = source_center, y = target_center, fill = value)
    ) +
      ggplot2::geom_tile() +
      ggplot2::facet_grid(
        rows = ggplot2::vars(projection_label),
        cols = ggplot2::vars(n),
        labeller = ggplot2::labeller(
          projection_label = ggplot2::label_value,
          n = ggplot2::label_both
        )
      ) +
      ggplot2::labs(
        title = title,
        x = "Source-coordinate center (u)",
        y = "Target-coordinate center (v)",
        fill = fill_label
      ) +
      ggplot2::coord_equal() +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    if (diverging) {
      limit <- max(abs(data$value), na.rm = TRUE)
      plot <- plot + ggplot2::scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0, limits = c(-limit, limit)
      )
    } else {
      plot <- plot + ggplot2::scale_fill_viridis_c(option = "C")
    }
    path <- file.path(plot_dir, paste0(prefix, filename))
    ggplot2::ggsave(path, plot, width = 12.5, height = 11, units = "in")
    path
  }

  list(
    representative = plot_heatmap(
      heatmap_data(representative_values, test_functions, n_values, "representative"),
      expression("Representative " * sqrt(n) *
                   integral(f * d * (hat(gamma)[n] - gamma[0]))),
      "Fluctuation", TRUE, "representative_fluctuations.pdf"
    ),
    mean = plot_heatmap(
      heatmap_data(mean_values, test_functions, n_values, "mean"),
      "Monte Carlo Mean of Scaled Fluctuations",
      "Mean", TRUE, "mean_fluctuations.pdf"
    ),
    sd = plot_heatmap(
      heatmap_data(sd_values, test_functions, n_values, "sd"),
      "Monte Carlo Standard Deviation of Scaled Fluctuations",
      "Standard deviation", FALSE, "sd_fluctuations.pdf"
    )
  )
}

select_diagnostic_features <- function(test_functions) {
  grid_size <- sqrt(nrow(test_functions$grid) /
                      nrow(test_functions$projection_spec))
  center_index <- ceiling(grid_size / 2)
  center_rows <- which(
    test_functions$grid$source_index == center_index &
      test_functions$grid$target_index == center_index
  )
  center_features <- test_functions$grid$feature[center_rows]
  center_columns <- match(center_features, test_functions$feature_names)
  names(center_columns) <- paste0(
    test_functions$grid$projection[center_rows], "_center"
  )
  c(total_mass = 1L, center_columns)
}

make_independent_reference <- function(checkpoint, population_integrals) {
  reference_n <- unique(checkpoint$diagnostics$n)
  if (length(reference_n) != 1L) {
    stop("The reference checkpoint must contain exactly one sample size.",
         call. = FALSE)
  }
  if (nrow(checkpoint$diagnostics) < 4L) {
    stop("At least four independent reference replications are needed.",
         call. = FALSE)
  }
  fluctuations <- compute_fluctuations(checkpoint, population_integrals)
  list(
    n = reference_n,
    replications = nrow(checkpoint$diagnostics),
    covariance = stats::cov(fluctuations),
    fluctuations = fluctuations
  )
}

make_gaussian_diagnostics <- function(checkpoint, fluctuations, test_functions,
                                      reference, plot_dir, prefix) {
  selected <- select_diagnostic_features(test_functions)
  selected_names <- names(selected)
  reference_sd <- sqrt(diag(reference$covariance))[selected]
  n_values <- sort(unique(checkpoint$diagnostics$n))

  histogram_pieces <- list()
  qq_pieces <- list()
  counter <- 0L
  for (n in n_values) {
    indices <- which(checkpoint$diagnostics$n == n)
    for (feature_index in seq_along(selected)) {
      column <- selected[[feature_index]]
      feature <- selected_names[[feature_index]]
      values <- fluctuations[indices, column]
      sd_reference <- reference_sd[[feature_index]]
      counter <- counter + 1L
      histogram_pieces[[counter]] <- data.frame(
        n = n, feature = feature, value = values,
        reference_sd = sd_reference
      )
      standardized <- sort(values / sd_reference)
      qq_pieces[[counter]] <- data.frame(
        n = n,
        feature = feature,
        theoretical = stats::qnorm(stats::ppoints(length(standardized))),
        observed = standardized
      )
    }
  }
  histogram_data <- do.call(rbind, histogram_pieces)
  qq_data <- do.call(rbind, qq_pieces)

  normal_overlay <- do.call(rbind, lapply(
    split(histogram_data, interaction(histogram_data$n, histogram_data$feature)),
    function(x) {
      limits <- range(x$value)
      grid <- seq(limits[[1]], limits[[2]], length.out = 200)
      data.frame(
        n = x$n[[1]], feature = x$feature[[1]], x = grid,
        density = stats::dnorm(grid, mean = 0, sd = x$reference_sd[[1]])
      )
    }
  ))

  histogram_plot <- ggplot2::ggplot(
    histogram_data, ggplot2::aes(x = value)
  ) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 25, fill = "grey80", color = "white"
    ) +
    ggplot2::geom_line(
      data = normal_overlay,
      ggplot2::aes(x = x, y = density),
      color = "#B2182B", linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    ggplot2::facet_grid(feature ~ n, scales = "free") +
    ggplot2::labs(
      title = paste0(
        "Scaled Fluctuations and Large-n Gaussian Reference (n = ",
        reference$n, ", independent R = ", reference$replications, ")"
      ),
      x = "Scaled fluctuation", y = "Density"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  qq_plot <- ggplot2::ggplot(
    qq_data, ggplot2::aes(x = theoretical, y = observed)
  ) +
    ggplot2::geom_point(size = 0.8, alpha = 0.65) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "#B2182B") +
    ggplot2::facet_grid(feature ~ n, scales = "free") +
    ggplot2::labs(
      title = "Normal Q--Q Plots Standardized by the Independent Reference Batch",
      x = "Standard normal quantile", y = "Observed standardized quantile"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  histogram_path <- file.path(plot_dir, paste0(prefix, "histograms.pdf"))
  qq_path <- file.path(plot_dir, paste0(prefix, "qqplots.pdf"))
  ggplot2::ggsave(histogram_path, histogram_plot, width = 12, height = 9, units = "in")
  ggplot2::ggsave(qq_path, qq_plot, width = 12, height = 9, units = "in")
  list(histograms = histogram_path, qqplots = qq_path)
}

make_covariance_diagnostics <- function(checkpoint, fluctuations, reference,
                                        plot_dir, prefix) {
  n_values <- sort(unique(checkpoint$diagnostics$n))
  reference_covariance <- reference$covariance
  reference_norm <- norm(reference_covariance, type = "F")
  reference_noise_benchmark <- sqrt(
    2 * (
      sum(reference_covariance^2) + sum(diag(reference_covariance))^2
    ) / (
      (reference$replications - 1) * sum(reference_covariance^2)
    )
  )
  comparison <- do.call(rbind, lapply(n_values, function(n) {
    indices <- which(checkpoint$diagnostics$n == n)
    covariance <- stats::cov(fluctuations[indices, , drop = FALSE])
    data.frame(
      n = n,
      relative_frobenius_error = norm(
        covariance - reference_covariance, type = "F"
      ) / reference_norm,
      reference_noise_benchmark = reference_noise_benchmark
    )
  }))

  stabilization_plot <- ggplot2::ggplot(
    comparison,
    ggplot2::aes(x = n, y = relative_frobenius_error)
  ) +
    ggplot2::geom_hline(
      yintercept = reference_noise_benchmark,
      linetype = "dashed", color = "#B2182B", linewidth = 0.7
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_x_continuous(breaks = n_values) +
    ggplot2::labs(
      title = "Empirical Covariance versus the Large-n Reference",
      subtitle = paste0(
        "Independent ", reference$replications,
        "-run reference at n = ", reference$n
      ),
      x = "Sample size per marginal (n)",
      y = "Relative Frobenius distance",
      caption = paste0(
        "Dashed line: Gaussian/Wishart plug-in same-law RMS benchmark = ",
        format(round(reference_noise_benchmark, 3), nsmall = 3)
      )
    ) +
    ggplot2::theme_minimal(base_size = 14)

  correlation <- stats::cov2cor(reference_covariance)
  correlation_data <- expand.grid(
    row = seq_len(nrow(correlation)), column = seq_len(ncol(correlation))
  )
  correlation_data$value <- as.vector(correlation)
  feature_names <- colnames(fluctuations)
  feature_groups <- ifelse(
    feature_names == "total_mass", "mass",
    sub("_bump.*$", "", feature_names)
  )
  group_runs <- rle(feature_groups)
  group_ends <- cumsum(group_runs$lengths)
  group_starts <- c(1, head(group_ends, -1) + 1)
  group_centers <- (group_starts + group_ends) / 2
  group_labels <- c(mass = "mass", p11 = "(1,1)", p12 = "(1,2)",
                    p21 = "(2,1)", p22 = "(2,2)")[group_runs$values]
  separators <- head(group_ends, -1) + 0.5
  correlation_plot <- ggplot2::ggplot(
    correlation_data,
    ggplot2::aes(x = column, y = row, fill = value)
  ) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1)
    ) +
    ggplot2::geom_vline(
      xintercept = separators, color = "grey25", linewidth = 0.25
    ) +
    ggplot2::geom_hline(
      yintercept = separators, color = "grey25", linewidth = 0.25
    ) +
    ggplot2::scale_x_continuous(
      breaks = group_centers, labels = group_labels,
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_reverse(
      breaks = group_centers, labels = group_labels,
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Reference Correlation of Scaled Fluctuations",
      subtitle = paste0(
        "Independent reference batch at n = ", reference$n,
        " (R = ", reference$replications, ")"
      ),
      x = "Projection block", y = "Projection block", fill = "Correlation"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  stabilization_path <- file.path(
    plot_dir, paste0(prefix, "covariance_stabilization.pdf")
  )
  correlation_path <- file.path(
    plot_dir, paste0(prefix, "reference_correlation.pdf")
  )
  ggplot2::ggsave(
    stabilization_path, stabilization_plot, width = 7.5, height = 5, units = "in"
  )
  ggplot2::ggsave(
    correlation_path, correlation_plot, width = 7, height = 6, units = "in"
  )
  list(
    stabilization = stabilization_path,
    correlation = correlation_path,
    comparison = comparison
  )
}

validate_population_against_empirical <- function(checkpoint, population) {
  largest_n <- max(checkpoint$diagnostics$n)
  observed <- checkpoint$diagnostics$plan_mass[
    checkpoint$diagnostics$n == largest_n
  ]
  relative_error <- abs(mean(observed) - population$mass) / population$mass
  if (length(observed) >= 20L && relative_error > 0.25) {
    warning(
      "At the largest n, mean empirical transported mass differs from the ",
      "closed-form population mass by ", round(100 * relative_error, 1),
      "%. Check solver and parameter conventions."
    )
  }
  invisible(relative_error)
}

main <- function() {
  args <- parse_arguments(commandArgs(trailingOnly = TRUE))
  assert_repository_root()
  require_packages(c("MASS", "ggplot2", "reticulate"))
  if (is.na(args$reps) || args$reps < 1L) {
    stop("--reps must be a positive integer.", call. = FALSE)
  }

  epsilon <- 0.5
  rho <- 1.0
  if (args$smoke_test) {
    n_grid <- c(20L, 30L)
    target_reps <- 4L
    grid_size <- 5L
    run_name <- "clt_gaussian_all_projections_smoke_test"
  } else {
    n_grid <- c(50L, 100L, 200L, 500L)
    target_reps <- args$reps
    grid_size <- 7L
    run_name <- "clt_gaussian_all_projections"
  }
  reference_n <- max(n_grid)
  target_reference_reps <- target_reps
  reference_seed_stream <- 1L

  marginals <- make_marginals()
  population <- gaussian_uot_population(
    marginals$a, marginals$b, marginals$A, marginals$B, epsilon, rho
  )
  test_functions <- make_test_functions(population, grid_size)

  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(args$plot_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- paste0(run_name, "_")
  checkpoint_path <- file.path(args$output_dir, paste0(prefix, "checkpoint.rds"))
  reference_checkpoint_path <- file.path(
    args$output_dir, paste0(prefix, "independent_reference_checkpoint.rds")
  )
  metadata_path <- file.path(args$output_dir, paste0(prefix, "metadata.rds"))
  reference_metadata_path <- file.path(
    args$output_dir, paste0(prefix, "independent_reference_metadata.rds")
  )
  summary_path <- file.path(args$output_dir, paste0(prefix, "summary.csv"))
  covariance_path <- file.path(
    args$output_dir, paste0(prefix, "independent_reference_covariance.rds")
  )

  metadata <- list(
    n_grid = n_grid,
    target_replications = target_reps,
    master_seed = args$seed,
    epsilon = epsilon,
    rho = rho,
    cost = "fixed squared Euclidean distance",
    entropy_reference = "KL(plan || empirical marginal product)",
    normalized_plan = FALSE,
    test_function_family = paste0(
      grid_size, "x", grid_size,
      " Gaussian bumps on all four source-target coordinate ",
      "projections (1,1), (1,2), (2,1), and (2,2), plus total mass"
    ),
    fluctuation_scaling = paste0(
      "sqrt(n), where n is the sample size of each marginal; ",
      "a total-sample convention differs by the constant sqrt(2)"
    ),
    marginals = marginals,
    population = population,
    test_functions = test_functions,
    created_at = Sys.time(),
    R_version = R.version.string,
    reticulate_python = Sys.getenv("RETICULATE_PYTHON", unset = NA_character_)
  )
  if (file.exists(metadata_path)) {
    previous <- readRDS(metadata_path)
    compatible <- identical(previous$n_grid, metadata$n_grid) &&
      identical(previous$master_seed, metadata$master_seed) &&
      identical(previous$epsilon, metadata$epsilon) &&
      identical(previous$rho, metadata$rho) &&
      identical(previous$marginals, metadata$marginals) &&
      identical(
        previous$test_functions$feature_names,
        metadata$test_functions$feature_names
      ) &&
      isTRUE(all.equal(
        previous$test_functions$population_integrals,
        metadata$test_functions$population_integrals,
        tolerance = 1e-12
      ))
    if (!compatible) {
      stop(
        "The existing run has a different grid, seed, parameters, marginals, ",
        "or test functions. Choose a different --output-dir before changing them.",
        call. = FALSE
      )
    }
    metadata$created_at <- previous$created_at
  }

  reference_metadata <- list(
    batch = "independent_large_n_reference",
    n = reference_n,
    target_replications = target_reference_reps,
    master_seed = args$seed,
    seed_stream = reference_seed_stream,
    seed_formula_version = "cell_seed_v2_stream_offset_1000000007",
    epsilon = epsilon,
    rho = rho,
    cost = "fixed squared Euclidean distance",
    entropy_reference = "KL(plan || empirical marginal product)",
    normalized_plan = FALSE,
    fluctuation_scaling = "sqrt(n), identical to the evaluation batch",
    marginals = marginals,
    population = population,
    feature_names = test_functions$feature_names,
    population_integrals = test_functions$population_integrals,
    created_at = Sys.time(),
    R_version = R.version.string,
    reticulate_python = Sys.getenv("RETICULATE_PYTHON", unset = NA_character_)
  )
  if (file.exists(reference_metadata_path)) {
    previous_reference <- readRDS(reference_metadata_path)
    compatible_reference <-
      identical(previous_reference$batch, reference_metadata$batch) &&
      identical(previous_reference$n, reference_metadata$n) &&
      identical(previous_reference$master_seed, reference_metadata$master_seed) &&
      identical(previous_reference$seed_stream, reference_metadata$seed_stream) &&
      identical(
        previous_reference$seed_formula_version,
        reference_metadata$seed_formula_version
      ) &&
      identical(previous_reference$epsilon, reference_metadata$epsilon) &&
      identical(previous_reference$rho, reference_metadata$rho) &&
      identical(previous_reference$marginals, reference_metadata$marginals) &&
      identical(previous_reference$feature_names, reference_metadata$feature_names) &&
      isTRUE(all.equal(
        previous_reference$population_integrals,
        reference_metadata$population_integrals,
        tolerance = 1e-12
      ))
    if (!compatible_reference) {
      stop(
        "The independent-reference metadata are incompatible with this run. ",
        "Choose a different --output-dir.", call. = FALSE
      )
    }
    reference_metadata$created_at <- previous_reference$created_at
  }

  checkpoint <- if (file.exists(checkpoint_path)) {
    readRDS(checkpoint_path)
  } else {
    empty_checkpoint(test_functions$feature_names)
  }
  reference_checkpoint <- if (file.exists(reference_checkpoint_path)) {
    readRDS(reference_checkpoint_path)
  } else {
    empty_checkpoint(test_functions$feature_names)
  }
  validate_checkpoint(
    checkpoint, test_functions$feature_names, n_grid,
    args$seed, 0L, "Evaluation"
  )
  validate_checkpoint(
    reference_checkpoint, test_functions$feature_names, reference_n,
    args$seed, reference_seed_stream, "Independent-reference"
  )
  atomic_save_rds(metadata, metadata_path)
  atomic_save_rds(reference_metadata, reference_metadata_path)
  full_design <- expand.grid(
    n = n_grid, replication = seq_len(target_reps), KEEP.OUT.ATTRS = FALSE
  )
  full_reference_design <- data.frame(
    n = rep(reference_n, target_reference_reps),
    replication = seq_len(target_reference_reps)
  )
  evaluation_seeds <- mapply(
    cell_seed, n = full_design$n, replication = full_design$replication,
    MoreArgs = list(master_seed = args$seed, seed_stream = 0L)
  )
  reference_seeds <- mapply(
    cell_seed,
    n = full_reference_design$n,
    replication = full_reference_design$replication,
    MoreArgs = list(
      master_seed = args$seed, seed_stream = reference_seed_stream
    )
  )
  if (length(intersect(evaluation_seeds, reference_seeds))) {
    stop("Evaluation and independent-reference seed streams overlap.",
         call. = FALSE)
  }

  design <- full_design
  if (nrow(checkpoint$diagnostics)) {
    completed <- paste(
      checkpoint$diagnostics$n, checkpoint$diagnostics$replication, sep = ":"
    )
    design_key <- paste(design$n, design$replication, sep = ":")
    design <- design[!design_key %in% completed, , drop = FALSE]
  }
  reference_design <- full_reference_design
  if (nrow(reference_checkpoint$diagnostics)) {
    reference_completed <- paste(
      reference_checkpoint$diagnostics$n,
      reference_checkpoint$diagnostics$replication,
      sep = ":"
    )
    reference_design_key <- paste(
      reference_design$n, reference_design$replication, sep = ":"
    )
    reference_design <- reference_design[
      !reference_design_key %in% reference_completed, , drop = FALSE
    ]
  }

  message(
    "Population transported mass: ", format(population$mass, digits = 6), "."
  )
  checkpoint <- complete_design(
    checkpoint, design, checkpoint_path, "Evaluation", 0L,
    args$seed, marginals, population, test_functions, epsilon, rho
  )
  reference_checkpoint <- complete_design(
    reference_checkpoint, reference_design, reference_checkpoint_path,
    "Independent reference", reference_seed_stream,
    args$seed, marginals, population, test_functions, epsilon, rho
  )
  validate_checkpoint(
    reference_checkpoint, test_functions$feature_names, reference_n,
    args$seed, reference_seed_stream, "Independent-reference"
  )

  selected <- checkpoint$diagnostics$n %in% n_grid &
    checkpoint$diagnostics$replication <= target_reps
  analysis_checkpoint <- list(
    diagnostics = checkpoint$diagnostics[selected, , drop = FALSE],
    integrals = checkpoint$integrals[selected, , drop = FALSE]
  )
  reference_selected <-
    reference_checkpoint$diagnostics$n == reference_n &
    reference_checkpoint$diagnostics$replication <= target_reference_reps
  analysis_reference_checkpoint <- list(
    diagnostics = reference_checkpoint$diagnostics[
      reference_selected, , drop = FALSE
    ],
    integrals = reference_checkpoint$integrals[
      reference_selected, , drop = FALSE
    ]
  )
  if (nrow(analysis_checkpoint$diagnostics) != nrow(full_design) ||
      nrow(analysis_reference_checkpoint$diagnostics) !=
        nrow(full_reference_design)) {
    stop("An evaluation or reference batch is incomplete.", call. = FALSE)
  }
  if (length(intersect(
    analysis_checkpoint$diagnostics$seed,
    analysis_reference_checkpoint$diagnostics$seed
  ))) {
    stop("Stored evaluation and reference seeds overlap.", call. = FALSE)
  }
  fluctuations <- compute_fluctuations(
    analysis_checkpoint, test_functions$population_integrals
  )
  validate_population_against_empirical(analysis_checkpoint, population)
  validate_population_against_empirical(
    analysis_reference_checkpoint, population
  )
  summary <- summarize_replications(
    analysis_checkpoint, fluctuations, population$mass
  )
  utils::write.csv(summary, summary_path, row.names = FALSE)

  if (target_reps >= 4L && target_reference_reps >= 4L) {
    reference <- make_independent_reference(
      analysis_reference_checkpoint, test_functions$population_integrals
    )
    atomic_save_rds(
      list(
        source = "independent_large_n_reference",
        n = reference$n,
        replications = reference$replications,
        master_seed = args$seed,
        seed_stream = reference_seed_stream,
        seeds = analysis_reference_checkpoint$diagnostics$seed,
        reference_replications =
          analysis_reference_checkpoint$diagnostics$replication,
        covariance = reference$covariance,
        feature_names = test_functions$feature_names
      ),
      covariance_path
    )
    heatmap_paths <- make_fluctuation_heatmaps(
      analysis_checkpoint, fluctuations, test_functions, args$plot_dir, prefix
    )
    gaussian_paths <- make_gaussian_diagnostics(
      analysis_checkpoint, fluctuations, test_functions, reference,
      args$plot_dir, prefix
    )
    covariance_paths <- make_covariance_diagnostics(
      analysis_checkpoint, fluctuations, reference, args$plot_dir, prefix
    )
    utils::write.csv(
      covariance_paths$comparison,
      file.path(args$output_dir, paste0(prefix, "covariance_stabilization.csv")),
      row.names = FALSE
    )
    invisible(c(heatmap_paths, gaussian_paths, covariance_paths[1:2]))
  }

  message("Saved checkpoint: ", checkpoint_path)
  message("Saved independent-reference checkpoint: ", reference_checkpoint_path)
  message("Saved metadata: ", metadata_path)
  message("Saved independent-reference metadata: ", reference_metadata_path)
  message("Saved summary: ", summary_path)
  if (target_reps >= 4L && target_reference_reps >= 4L) {
    message("Saved diagnostic plots in: ", args$plot_dir)
  }
}

if (sys.nframe() == 0L) main()
