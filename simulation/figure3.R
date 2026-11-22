# Run from the repository root with Rscript --vanilla simulation/figure3.R
source("R/uniform_kl_log_plan.R", local = TRUE)
source("R/thresh_ot_func.R", local = TRUE)

# The published example uses normalized costs in [0,1] and epsilon = 0.01,
# so the Gibbs kernel is representable in double precision.
figure3_balanced_plan <- function(cost, epsilon) {
  kernel <- exp(-cost / epsilon)
  a <- rep(1 / nrow(cost), nrow(cost))
  b <- rep(1 / ncol(cost), ncol(cost))
  v <- rep(1, length(b))
  for (iteration in seq_len(100000L)) {
    u <- a / as.vector(kernel %*% v)
    v <- b / as.vector(crossprod(kernel, u))
    if (iteration %% 100L == 0L) {
      plan <- kernel * outer(u, v)
      error <- max(abs(rowSums(plan) - a), abs(colSums(plan) - b))
      if (!is.finite(error)) stop("Nonfinite balanced marginal residual.")
      if (error <= 1e-12) return(plan)
    }
  }
  stop("Balanced Sinkhorn did not reach the marginal tolerance.")
}

fit_figure3 <- function() {
  sample <- read.csv("simulation/ET_round2/figure3/sample.csv")
  sample <- sample[order(sample$id), ]
  control <- which(sample$group == "control")
  treated <- which(sample$group == "treated")
  x <- as.matrix(sample[c("x1", "x2")])
  stopifnot(identical(sample$id, seq_len(200L)), length(control) == 100L,
            length(treated) == 100L, all(is.finite(x)))
  # Plan rows correspond to treated observations and columns to controls.
  distances <- as.matrix(dist(x))[treated, control]
  cost <- distances / max(distances)
  epsilon <- 0.01
  rho <- c(0.002, 0.01, Inf)
  plans <- lapply(rho[1:2], function(r) {
    exp(uniform_kl_log_plan(cost, epsilon, r)$log_plan)
  })
  plans[[3]] <- figure3_balanced_plan(cost, epsilon)
  list(x = x, control = control, treated = treated, plans = plans,
       epsilon = epsilon, rho = rho)
}

draw_figure3 <- function(fit) {
  previous <- par(mfrow = c(1, 3))
  on.exit(par(previous))
  titles <- c(paste("UOT: rho=", fit$rho[1:2]), "OT")
  panels <- lapply(seq_along(fit$plans), function(i) {
    plan <- fit$plans[[i]]
    pairs <- add_top_match(plan, fit$x, fit$control, fit$treated, 1:3,
                           main = titles[i], cex = 1.2)
    data.frame(epsilon = fit$epsilon, rho = fit$rho[i], rank = 1:3,
               treated_index = pairs[, 1], control_index = pairs[, 2],
               treated_id = fit$treated[pairs[, 1]], control_id = fit$control[pairs[, 2]],
               matched_mass = plan[pairs])
  })
  do.call(rbind, panels)
}

write_figure3 <- function(out) {
  fit <- fit_figure3()
  if (!dir.exists(out) && !dir.create(out, recursive = TRUE)) stop("Cannot create ", out)
  pdf(file.path(out, "rho_top3.pdf"), width = 7, height = 4, useDingbats = FALSE)
  ranks <- tryCatch(draw_figure3(fit), finally = dev.off())
  write.csv(ranks, file.path(out, "ranked_pairs.csv"), row.names = FALSE)
  invisible(ranks)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L) stop("Usage: Rscript simulation/figure3.R [output-directory]")
  write_figure3(if (length(args)) args[[1L]] else "plot/figure3")
}
