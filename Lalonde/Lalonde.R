# Reproduce Table 1 from the converged Figure 5 grid and the original cohort.
# Run Rscript Lalonde/Lalonde.R --uot-only for the three primary contrasts.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) && !identical(args, "--uot-only"))
  stop("Usage: Rscript Lalonde/Lalonde.R [--uot-only]")
uot_only <- identical(args, "--uot-only")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
# Some Rscript launchers encode spaces in the --file argument.
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
folder <- dirname(normalizePath(script, mustWork = TRUE))
if (!uot_only) {
  packages <- c("MatchIt", "CBPS", "sbw", "quadprog")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install the Table 1 comparison packages: ", paste(missing, collapse = ", "))
  if (packageVersion("CBPS") != "0.23")
    stop("Table 1 uses CBPS 0.23. See Lalonde/README.md for installation instructions.")
}

# One run supplies the same certified grid, data and selected cell to both figures and Table 1.
status <- system2(file.path(R.home("bin"), "Rscript"),
                  c("--vanilla", shQuote(file.path(folder, "UOT_sensi.R"))))
if (status != 0L) stop("Could not reproduce the converged LaLonde grid.")
output <- file.path(folder, "results")
grid <- read.csv(file.path(output, "converged_grid.csv"))
selected <- subset(grid, is_selected_table1 == 1)
stopifnot(nrow(selected) == 1L, nrow(grid) == 300L, all(grid$max_abs_foc < 1e-10))
message(sprintf("Selected Figure 5 cell: epsilon=%g, rho=%g", selected$epsilon, selected$rho))

covars <- c("age", "educ", "black", "hisp", "married", "nodegree", "re75")
columns <- c("treat", covars, "re78")
read_sample <- function(name, psid = FALSE) {
  dat <- read.table(file.path(folder, "data", paste0(name, ".txt")))
  if (psid) dat <- dat[, -8L]  # Only the three PSID files contain re74.
  setNames(dat, columns)
}
treated <- read_sample("nsw_treated")
randomized <- read_sample("nsw_control")
controls <- do.call(rbind, lapply(c("psid_controls", "psid2_controls", "psid3_controls"),
                                  read_sample, psid = TRUE))
observed <- rbind(treated, controls)
unadjusted <- mean(treated$re78) - mean(controls$re78)
experimental <- mean(treated$re78) - mean(randomized$re78)
stopifnot(abs(selected$experimental_benchmark - experimental) < 1e-8)

estimates <- c(Unadjusted = unadjusted)
if (!uot_only) {
  source(file.path(folder, "lalonde_help.R"))
  model <- reformulate(covars, "treat")
  matched_contrast <- function(fit) {
    with(MatchIt::match_data(fit, data = observed),
         mean(re78[treat == 1]) - mean(re78[treat == 0]))
  }
  weighted_contrast <- function(w) {
    mean(treated$re78) - weighted.mean(controls$re78, w[observed$treat == 0])
  }
  # Appendix B specifies a standardized 0.1 caliper for the headline row.
  nearest <- MatchIt::matchit(model, data = observed, method = "nearest",
                              distance = "glm", caliper = 0.1, std.caliper = TRUE)
  cbps <- CBPS::CBPS(model, data = observed, ATT = 1)
  maha <- MatchIt::matchit(model, data = observed, method = "nearest", distance = "mahalanobis")
  stable <- sbw_att(observed, "treat", "re78", covars, bal_std = "target",
                    bal_tol = 0.02, bal_alg = FALSE)
  # The published overlap comparison uses the raw treated mean and
  # e*(1-e)-weighted controls, with e fitted by logistic regression.
  ps <- fitted(glm(model, data = observed, family = binomial()))
  estimates <- c(estimates, "PS Nearest" = matched_contrast(nearest),
                 CBPS = weighted_contrast(cbps$weights), Mahalanobis = matched_contrast(maha),
                 SBW = stable$att, "Overlap Weights" = weighted_contrast(ps * (1 - ps)))
}
estimates <- c(estimates, UOT = selected$estimate, Experimental = experimental)
stopifnot(all(is.finite(estimates)))
table1 <- data.frame(method = names(estimates), estimate = unname(estimates))
filename <- if (uot_only) "table1_uot.csv" else "table1.csv"
write.csv(table1, file.path(output, filename), row.names = FALSE)
print(transform(table1, estimate = round(estimate, 1)), row.names = FALSE)
