# Product-reference KL transport between equally weighted observations.
# Return the raw log plan even when individual masses underflow in ordinary arithmetic.
uniform_kl_log_plan <- function(cost, epsilon, rho, max_iter = 2000L) {
  stopifnot(is.matrix(cost), nrow(cost) > 0L, ncol(cost) > 0L,
            all(is.finite(cost)), all(cost >= 0),
            length(epsilon) == 1L, is.finite(epsilon), epsilon > 0,
            length(rho) == 1L, is.finite(rho), rho > 0)
  log_a <- -log(nrow(cost))
  log_b <- -log(ncol(cost))
  tau <- rho / (rho + epsilon)
  log_kernel <- -cost / epsilon
  row_lse <- function(x) {
    peak <- apply(x, 1L, max)
    peak + log(rowSums(exp(x - peak)))
  }

  log_v <- numeric(ncol(cost))
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
    if (!is.finite(log_mass) || log_mass > 1e-10) stop("Invalid raw UOT mass.")
    return(list(log_plan = log_plan, log_mass = log_mass,
                max_abs_foc = residual, iterations = iteration))
  }
  stop("Log-domain Sinkhorn did not meet the first-order tolerance.")
}
