# Product-reference KL transport between two probability histograms. The
# ground cost and both marginal penalties are in the manuscript's units.
acs_kl_plan <- function(a, b, cost, eps, rho, tol = 1e-10) {
  stopifnot(is.numeric(a), is.numeric(b), is.numeric(cost), is.matrix(cost),
            nrow(cost) == length(a), ncol(cost) == length(b),
            all(is.finite(cost)), all(cost >= 0),
            all(is.finite(a)), all(a >= 0), abs(sum(a) - 1) < 1e-12,
            all(is.finite(b)), all(b >= 0), abs(sum(b) - 1) < 1e-12,
            length(eps) == 1L, is.finite(eps), eps > 0,
            length(rho) == 1L, is.finite(rho), rho > 0,
            length(tol) == 1L, is.finite(tol), tol > 0)

  # Entries outside the support of a x b are exactly zero at finite KL.
  active_a <- a > 0
  active_b <- b > 0
  ground <- cost[active_a, active_b, drop = FALSE]
  log_a <- log(a[active_a])
  log_b <- log(b[active_b])
  tau <- rho / (rho + eps)
  row_logsum <- function(x) {
    peak <- apply(x, 1L, max)
    peak + log(rowSums(exp(x - peak)))
  }

  # Given source potential f, the target stationarity equation determines g.
  # The remaining log-domain equation is f + tau*eps*row_logsum(scores) = 0.
  evaluate <- function(f) {
    from_source <- (f - ground) / eps + log_a
    peaks <- apply(from_source, 2L, max)
    centered <- sweep(from_source, 2L, peaks, "-")
    log_normalizers <- log(colSums(exp(centered)))
    log_given_target <- sweep(centered, 2L, log_normalizers, "-")
    col_sums <- peaks + log_normalizers
    g <- -tau * eps * col_sums
    from_target <- sweep(-ground / eps, 2L, g / eps + log_b, "+")
    row_sums <- row_logsum(from_target)
    # Store each column using its conditional probabilities and log mass.
    # This avoids cancellation between f and g at large rho/eps.
    log_col_mass <- log_b + eps / (rho + eps) * col_sums
    list(error = f + tau * eps * row_sums,
         log_plan = sweep(log_given_target, 2L, log_col_mass, "+"),
         given_source = exp(from_target - row_sums),
         given_target = exp(log_given_target))
  }

  f <- numeric(length(log_a))
  state <- evaluate(f)
  for (iteration in seq_len(200L)) {
    log_plan <- state$log_plan
    log_ratio <- log_plan - outer(log_a, log_b, "+")
    row_penalty <- rho * (row_logsum(log_plan) - log_a)
    col_penalty <- rho * (row_logsum(t(log_plan)) - log_b)
    # Check the full cellwise derivative of the original primal objective.
    residual <- max(abs(ground + eps * log_ratio +
                          outer(row_penalty, col_penalty, "+")))
    if (is.finite(residual) && residual <= tol) {
      plan <- matrix(0, nrow(cost), ncol(cost), dimnames = dimnames(cost))
      plan[active_a, active_b] <- exp(log_plan)
      mass <- sum(plan)
      if (!is.finite(mass) || mass <= 0 || mass > 1 + tol) {
        stop("Certified ACS plan violates the total-mass bound: ", mass)
      }
      attr(plan, "max_abs_foc") <- residual
      attr(plan, "transport_mass") <- mass
      attr(plan, "iterations") <- iteration
      return(plan)
    }

    # The derivative after eliminating g is I - tau^2 P(j|i) P(i|j).
    jacobian <- diag(length(f)) -
      tau^2 * tcrossprod(state$given_source, state$given_target)
    step <- drop(solve(jacobian, state$error))
    rate <- 1
    repeat {
      proposal <- f - rate * step
      candidate <- evaluate(proposal)
      if (max(abs(candidate$error)) <=
          (1 - 1e-4 * rate) * max(abs(state$error))) break
      rate <- rate / 2
      if (rate < 2^-20) {
        # One alternating log-Sinkhorn step is always a contraction.
        proposal <- f - state$error
        candidate <- evaluate(proposal)
        break
      }
    }
    f <- proposal
    state <- candidate
  }
  stop(sprintf("ACS KL solve did not converge (eps=%g, rho=%g, residual=%.3g)",
               eps, rho, residual))
}
