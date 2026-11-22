# Jarque–Bera diagnostics on the same evaluation IDs at each sample size.
jb_sample_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)
jb_features <- c(
  "total_mass", "p11_bump_s04_t04", "p12_bump_s04_t04",
  "p21_bump_s04_t04", "p22_bump_s04_t04"
)

compute_jarque_bera_table <- function(diagnostics, fluctuations, count = 3000L) {
  if (length(count) != 1L || !is.finite(count) || count != as.integer(count) ||
      count < 4L || count > 3000L) {
    stop("The JB replication count must be an integer from 4 to 3000.")
  }
  if (!all(c("n", "replication") %in% names(diagnostics)) ||
      nrow(fluctuations) != nrow(diagnostics) ||
      !all(jb_features %in% colnames(fluctuations))) {
    stop("Misaligned JB input or missing diagnostic features.")
  }
  ids <- seq_len(count)
  table <- do.call(rbind, lapply(jb_sample_sizes, function(n) {
    rows <- which(diagnostics$n == n & diagnostics$replication %in% ids)
    indices <- match(ids, diagnostics$replication[rows])
    if (length(rows) != count || anyNA(indices)) {
      stop(sprintf("Expected exactly one copy of evaluation IDs 1:%d at n=%d.", count, n))
    }
    values <- fluctuations[rows[indices], jb_features, drop = FALSE]
    if (any(!is.finite(values))) stop("Nonfinite fluctuation in JB input.")
    centered <- sweep(values, 2L, colMeans(values), "-")
    second <- colMeans(centered^2)
    if (any(second == 0)) stop("JB requires positive variance for every feature.")
    skewness <- colMeans(centered^3) / second^(3/2)
    excess <- colMeans(centered^4) / second^2 - 3
    jb <- count / 6 * (skewness^2 + excess^2 / 4)
    data.frame(
      n = n, feature = jb_features, replications = count, jb = unname(jb),
      p_value = unname(stats::pchisq(jb, df = 2, lower.tail = FALSE))
    )
  }))
  table$p_holm <- stats::p.adjust(table$p_value, method = "holm")
  table
}

verify_published_jarque_bera <- function(table) {
  # Entries follow increasing n and the feature order above. Preserve trailing
  # zeros in the p-values because each cell has its own published precision.
  jb_target <- c(
    121.0, 171.3, 301.9, 324.8, 168.6,
    116.5,  58.6, 125.7,  84.9, 118.8,
     60.3,  16.2,  48.0,  29.8,  22.9,
     26.7,   4.1,  14.3,  22.6,  20.8,
     10.0,   2.4,   8.4,   6.1,   6.0,
      5.1,   8.6,  10.9,   8.5,   4.5
  )
  holm_target <- c(
    "1.4e-25", "1.7e-36", "8.1e-65", "8.9e-70", "6.6e-36",
    "1.2e-24", "3.7e-12", "1.3e-26", "8.0e-18", "3.8e-25",
    "1.7e-12", "0.0039",  "7.2e-10", "6.1e-6",  "0.0002",
    "2.7e-5",  "0.320",   "0.0095",  "0.0002",  "0.0004",
    "0.0688",  "0.320",   "0.120",   "0.283",   "0.283",
    "0.320",   "0.120",   "0.0474",  "0.120",   "0.320"
  )
  if (!identical(as.integer(table$n), rep(jb_sample_sizes, each = length(jb_features))) ||
      !identical(as.character(table$feature), rep(jb_features, length(jb_sample_sizes))) ||
      any(table$replications != 3000L)) {
    stop("The published JB comparison requires 3000 replications in all 30 cells.")
  }
  decimals <- nchar(sub(".*[.]", "", sub("e.*", "", holm_target)))
  formats <- ifelse(grepl("e", holm_target), "%.*e", "%.*f")
  rounded_p <- as.numeric(sprintf(formats, decimals, table$p_holm))
  bad <- which(!is.finite(table$jb) | !is.finite(table$p_holm) |
                 round(table$jb, 1) != jb_target | rounded_p != as.numeric(holm_target))
  if (length(bad)) {
    stop("Published JB / Holm comparison failed:\n", paste(sprintf(
      "n=%d %s: computed %.17g / %.17g; expected %.1f / %s",
      table$n[bad], table$feature[bad], table$jb[bad], table$p_holm[bad],
      jb_target[bad], holm_target[bad]
    ), collapse = "\n"), call. = FALSE)
  }
  message("All 30 JB statistics and all 30 Holm p-values match at the published precision.")
  invisible(TRUE)
}

write_jarque_bera_table <- function(diagnostics, fluctuations, output_dir,
                                    count = 3000L, verify = FALSE) {
  table <- compute_jarque_bera_table(diagnostics, fluctuations, count)
  if (verify) verify_published_jarque_bera(table)
  path <- file.path(output_dir, sprintf("jarque_bera_first%d.csv", count))
  # Seventeen significant digits round-trip IEEE double values through CSV.
  output <- table
  for (field in c("jb", "p_value", "p_holm")) {
    output[[field]] <- sprintf("%.17g", table[[field]])
  }
  utils::write.table(output, path, sep = ",", row.names = FALSE, quote = FALSE)
  invisible(path)
}
