#!/usr/bin/env Rscript
# Run from the repository root: Rscript --vanilla simulation/MR_R1Q2_2.R
source("R/uniform_kl_log_plan.R", local = TRUE)
source("simulation/MR_functions.R", local = TRUE)

# Strict inequality implements "farther than the radius". Evaluate the raw
# mass in log space so small fitted entries are included in the sum.
log_mass_beyond <- function(log_plan, distance, radius) {
  far <- log_plan[distance > radius]
  if (!length(far)) return(-Inf)
  peak <- max(far)
  peak + log(sum(exp(far - peak)))
}

simulate_figure2 <- function() {
  params <- get_marg_norm(2)
  n <- 100L
  radius <- 0.1
  penalties <- c(0.1, 0.05, 0.02, 0.01, 0.005, 0.002)
  draws <- lapply(seq_len(50L), function(replication) {
    seed <- 20261121L + replication
    set.seed(seed)
    x <- MASS::mvrnorm(n, params$a, params$A)
    y <- MASS::mvrnorm(n, params$b, params$B)
    squared <- outer(rowSums(x*x), rowSums(y*y), "+") - 2 * tcrossprod(x, y)
    distance <- sqrt(pmax(squared, 0))
    settings <- lapply(penalties, function(t) {
      fit <- uniform_kl_log_plan(distance, epsilon = t, rho = t)
      log_far <- log_mass_beyond(fit$log_plan, distance, radius)
      log_bound <- -radius / (t + 2*t)
      if (log_far > min(log_bound, fit$log_mass) + 1e-10)
        stop(sprintf("Off-radius mass check failed at seed %d, epsilon=rho=%g", seed, t))
      data.frame(replication, seed, epsilon = t, rho = t, radius,
                 off_radius_mass = exp(log_far), total_mass = exp(fit$log_mass),
                 upper_bound = exp(log_bound), log_off_radius_mass = log_far,
                 log_total_mass = fit$log_mass, log_upper_bound = log_bound,
                 max_abs_foc = fit$max_abs_foc, iterations = fit$iterations)
    })
    do.call(rbind, settings)
  })
  do.call(rbind, draws)
}

summarize_figure2 <- function(draws) {
  settings <- lapply(unique(draws$epsilon), function(t) {
    rows <- draws[draws$epsilon == t, ]
    q <- unname(quantile(rows$off_radius_mass, c(0.25, 0.5, 0.75)))
    data.frame(epsilon = t, rho = rows$rho[1], radius = rows$radius[1],
               replications = nrow(rows), q25 = q[1], median = q[2], q75 = q[3],
               upper_bound = rows$upper_bound[1], median_total_mass = median(rows$total_mass))
  })
  do.call(rbind, settings)
}

plot_figure2 <- function(summary) {
  library(ggplot2)
  ggplot(summary, aes(epsilon, group = 1)) +
    geom_line(aes(y = upper_bound, linetype = "Upper bound"), linewidth = 0.6) +
    geom_line(aes(y = median, linetype = "Median"), linewidth = 0.6) +
    geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.045, linewidth = 0.5) +
    geom_point(aes(y = median), size = 1.7) +
    scale_x_continuous(transform = scales::transform_compose("log10", "reverse"),
                       breaks = summary$epsilon, labels = as.character(summary$epsilon)) +
    scale_linetype_manual(values = c("Median" = "solid", "Upper bound" = "dashed"),
                          breaks = c("Median", "Upper bound"), name = NULL) +
    expand_limits(y = c(1e-10, 1)) +
    scale_y_log10(breaks = 10^seq(0, -10, by = -2),
                  labels = scales::label_log()) +
    labs(x = expression("Joint penalty  " * epsilon == rho),
         y = "Matched mass beyond distance 0.1") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L) stop("Usage: Rscript simulation/MR_R1Q2_2.R [output-directory]")
  out <- if (length(args)) args[[1L]] else "simulation/ET_round2/figure2"
  draws <- simulate_figure2()
  stats <- summarize_figure2(draws)
  stopifnot(nrow(draws) == 300L, nrow(stats) == 6L, all(stats$replications == 50L))
  if (!dir.exists(out) && !dir.create(out, recursive = TRUE)) stop("Cannot create ", out)
  write.csv(draws, file.path(out, "replications.csv"), row.names = FALSE)
  write.csv(stats, file.path(out, "summary.csv"), row.names = FALSE)
  ggplot2::ggsave(file.path(out, "rho.pdf"), plot_figure2(stats), width = 6.8,
                  height = 3.6, units = "in", device = grDevices::pdf, useDingbats = FALSE)
  message(sprintf("Saved results for %d fits. Maximum absolute residual %.3g.",
                  nrow(draws), max(draws$max_abs_foc)))
}
