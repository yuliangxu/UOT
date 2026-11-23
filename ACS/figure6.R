# Run from the repository root: Rscript --vanilla ACS/figure6.R [input.csv.gz] [output-directory]
source("ACS/acs_kl_solver.R")

# Extract the paired ages needed for Figure 6, in the couple-ID order used by
# the original dplyr bootstrap. No individual records are written to disk.
acs_couple_ages <- function(path) {
  if (!file.exists(path)) stop("ACS input file does not exist: ", path)
  required <- c("YEAR", "SERIAL", "PERNUM", "SPLOC", "SEX", "AGE")
  fields <- names(utils::read.csv(path, nrows = 0L))
  missing <- setdiff(required, fields)
  if (length(missing)) stop("Missing ACS columns: ", paste(missing, collapse = ", "))
  classes <- ifelse(fields %in% required, "integer", "NULL")
  people <- utils::read.csv(path, colClasses = classes)
  if (length(unique(people$YEAR)) != 1L || anyNA(people$YEAR))
    stop("Figure 6 requires a single ACS year because household serial numbers repeat across years.")
  year <- people$YEAR[1L]
  people <- people[!is.na(people$SPLOC) & people$SPLOC > 0 & people$SEX %in% 1:2, ]
  if (anyNA(people[c("SERIAL", "PERNUM")])) stop("Missing household or person identifier.")
  id <- with(people, paste(SERIAL, pmin(PERNUM, SPLOC), pmax(PERNUM, SPLOC), sep = "_"))
  # Keep linked pairs with two records, one for each spouse.
  ord <- order(id, method = "radix")
  groups <- rle(id[ord])
  ends <- cumsum(groups$lengths)
  last <- ends[groups$lengths == 2L]
  first_row <- ord[last - 1L]
  last_row <- ord[last]
  opposite <- people$SEX[first_row] != people$SEX[last_row]
  first_row <- first_row[opposite]
  last_row <- last_row[opposite]
  first_is_male <- people$SEX[first_row] == 1L
  male <- people$AGE[ifelse(first_is_male, first_row, last_row)]
  female <- people$AGE[ifelse(first_is_male, last_row, first_row)]
  if (!length(male) || anyNA(c(male, female))) stop("No complete paired age cohort.")
  list(male = male, female = female, year = year)
}

figure6_bootstrap <- function(pairs) {
  k <- 30L
  age_range <- range(pairs$male, pairs$female)
  breaks <- seq(age_range[1], age_range[2], length.out = k + 1L)
  centers <- (breaks[-1L] + breaks[-length(breaks)]) / 2
  cost <- abs(outer(centers, centers, "-"))
  cost <- cost / max(cost)
  male_bins <- cut(pairs$male, breaks, include.lowest = TRUE, labels = FALSE)
  female_bins <- cut(pairs$female, breaks, include.lowest = TRUE, labels = FALSE)
  seeds <- 124:173
  n_pairs <- 1000L
  margins <- lapply(seeds, function(seed) {
    set.seed(seed)
    draw <- sample.int(length(male_bins), n_pairs, replace = TRUE)
    list(male = tabulate(male_bins[draw], k) / n_pairs,
         female = tabulate(female_bins[draw], k) / n_pairs)
  })
  grid <- expand.grid(eps = c(0.001, 0.01, 0.1), rho = c(0.01, 0.1, 1, 10))
  cells <- expand.grid(male_bin = seq_len(k), female_bin = seq_len(k))
  panels <- diagnostics <- vector("list", nrow(grid))
  for (g in seq_len(nrow(grid))) {
    eps <- grid$eps[g]
    rho <- grid$rho[g]
    fits <- lapply(margins, function(h) acs_kl_plan(h$male, h$female, cost, eps, rho))
    diagnostics[[g]] <- data.frame(eps, rho, seed = seeds,
      male_bins_used = vapply(margins, function(h) sum(h$male > 0), integer(1)),
      female_bins_used = vapply(margins, function(h) sum(h$female > 0), integer(1)),
      iterations = vapply(fits, function(fit) as.numeric(attr(fit, "iterations")), numeric(1)),
      max_abs_foc = vapply(fits, function(fit) attr(fit, "max_abs_foc"), numeric(1)),
      total_mass = vapply(fits, sum, numeric(1)))
    plans <- vapply(fits, as.vector, numeric(k * k))
    quantiles <- apply(plans, 1L, stats::quantile, probs = c(.025, .975), names = FALSE)
    panels[[g]] <- data.frame(eps, rho, cells, mean = rowMeans(plans),
                              q025 = quantiles[1L, ], q975 = quantiles[2L, ])
  }
  histograms <- do.call(rbind, Map(function(seed, h) {
    data.frame(seed, bin = seq_len(k), lower = breaks[-length(breaks)], upper = breaks[-1L],
               male = h$male, female = h$female)
  }, seeds, margins))
  list(panels = do.call(rbind, panels), diagnostics = do.call(rbind, diagnostics),
       histograms = histograms, breaks = breaks)
}

figure6_plot <- function(result) {
  library(ggplot2)
  p <- result$panels
  stat_names <- c("Mean", "2.5%", "97.5%")
  tile <- p[rep(seq_len(nrow(p)), times = 3L), c("eps", "rho", "male_bin", "female_bin")]
  tile$stat <- factor(rep(stat_names, each = nrow(p)), levels = stat_names)
  tile$mass <- c(p$mean, p$q025, p$q975)
  breaks <- result$breaks
  labs <- paste0(breaks[-length(breaks)], "–", breaks[-1L])
  tile$male_bin <- factor(tile$male_bin, seq_along(labs), labels = labs)
  tile$female_bin <- factor(tile$female_bin, seq_along(labs), labels = labs)
  ticks <- labs[seq(1, length(labs), by = 5)]
  pow035 <- scales::new_transform("pow035", function(x) x ^ .35, function(x) x ^ (1/.35))
  ggplot(tile, aes(female_bin, male_bin, fill = mass)) +
    geom_tile() +
    scale_x_discrete(breaks = ticks) + scale_y_discrete(breaks = ticks) +
    facet_grid(rho ~ eps + stat, labeller = labeller(
      rho = function(x) paste0("rho = ", x), eps = function(x) paste0("eps = ", x))) +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, max(tile$mass)),
                        transform = pow035, name = NULL) +
    ggplot2::labs(title = "UOT bootstrap summaries (var=AGE, n_pairs=1000, B=50)",
                  x = "Female age bin", y = "Male age bin") +
    theme_minimal() + theme(strip.text = element_text(size = 9),
                            axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                            axis.text.y = element_text(size = 7))
}

run_figure6 <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) > 2L) stop("Usage: Rscript ACS/figure6.R [input.csv.gz] [output-directory]")
  input <- if (length(args)) args[1L] else "data/usa_00002.csv.gz"
  out <- if (length(args) == 2L) args[2L] else "ACS/results"
  pairs <- acs_couple_ages(input)
  message(sprintf("ACS %d: %d couples. Fitting 50 replicates at each of 12 settings.", pairs$year, length(pairs$male)))
  result <- figure6_bootstrap(pairs)
  if (!dir.exists(out) && !dir.create(out, recursive = TRUE)) stop("Cannot create output directory: ", out)
  ggplot2::ggsave(file.path(out, "acs.png"), figure6_plot(result), device = grDevices::png,
                  width = 1399, height = 588, units = "px", dpi = 100, bg = "white")
  for (name in c("panels", "diagnostics", "histograms"))
    utils::write.csv(result[[name]], file.path(out, paste0("figure6_", name, ".csv")), row.names = FALSE)
  message(sprintf("Wrote Figure 6 and aggregate checks. Maximum absolute FOC residual %.3g.",
                  max(result$diagnostics$max_abs_foc)))
}

if (sys.nframe() == 0L) run_figure6()
