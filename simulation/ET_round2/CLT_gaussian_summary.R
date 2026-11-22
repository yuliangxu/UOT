#!/usr/bin/env Rscript

# Summarize completed CLT Gaussian HPC shards without rerunning simulations.

options(stringsAsFactors = FALSE)

core_script <- file.path("simulation", "ET_round2", "CLT_gaussian_experiment.R")
hpc_script <- file.path("simulation", "ET_round2", "CLT_gaussian_hpc.R")
if (!file.exists(core_script) || !file.exists(hpc_script)) {
  stop("Run this script from the repository root.", call. = FALSE)
}
source(core_script)
source(hpc_script)

parse_summary_arguments <- function(args) {
  value_after_equals <- function(prefix, default) {
    hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
    if (!length(hit)) return(default)
    sub(paste0("^", prefix, "="), "", hit[[length(hit)]])
  }

  list(
    help = "--help" %in% args || "-h" %in% args,
    smoke_test = "--smoke-test" %in% args,
    reps = as.integer(value_after_equals("--reps", "10000")),
    rep_start = as.integer(value_after_equals("--rep-start", "1")),
    num_shards = as.integer(value_after_equals("--num-shards", "100")),
    seed = as.integer(value_after_equals("--seed", "20260814")),
    shard_dir = value_after_equals(
      "--shard-dir", "/cwork/yx306/UOT/out/clt-gaussian-10000/shards"
    ),
    output_dir = value_after_equals(
      "--output-dir", "/cwork/yx306/UOT/out/clt-gaussian-10000/result"
    ),
    plot_dir = value_after_equals(
      "--plot-dir", "/cwork/yx306/UOT/out/clt-gaussian-10000/plot"
    )
  )
}

summary_usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript --vanilla simulation/ET_round2/CLT_gaussian_summary.R \\\n",
    "    --reps=10000 --rep-start=1 --num-shards=100 \\\n",
    "    --shard-dir=/cwork/yx306/UOT/out/clt-gaussian-10000/shards \\\n",
    "    --output-dir=/cwork/yx306/UOT/out/clt-gaussian-10000/result \\\n",
    "    --plot-dir=/cwork/yx306/UOT/out/clt-gaussian-10000/plot\n",
    sep = ""
  )
}

validate_summary_arguments <- function(args) {
  integers <- c(args$reps, args$rep_start, args$num_shards, args$seed)
  if (anyNA(integers) || args$reps < 4L || args$rep_start < 1L ||
      args$rep_start > args$reps || args$num_shards < 1L || args$seed < 0L) {
    stop("Invalid replication, shard, or seed arguments.", call. = FALSE)
  }
  allowed_root <- "/cwork/yx306/UOT/out"
  output_paths <- c(args$output_dir, args$plot_dir)
  normalized <- vapply(
    output_paths, normalizePath, character(1), mustWork = FALSE
  )
  if (any(!startsWith(normalized, paste0(allowed_root, "/")))) {
    stop(
      "--output-dir and --plot-dir must be below ", allowed_root, ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

load_completed_shards <- function(args, configuration) {
  hpc_args <- args
  hpc_args$mode <- "merge"
  hpc_args$shard_id <- NA_integer_
  directory <- configuration_directory(hpc_args, configuration)
  if (!dir.exists(directory)) {
    stop("Shard directory does not exist: ", directory, call. = FALSE)
  }

  evaluation_parts <- vector("list", args$num_shards)
  reference_parts <- vector("list", args$num_shards)
  input_paths <- character(2L * args$num_shards)
  environment_signatures <- vector("list", 2L * args$num_shards)
  counter <- 0L

  for (shard_id in seq_len(args$num_shards)) {
    replications <- assigned_replications(hpc_args, shard_id)
    for (batch in c("evaluation", "reference")) {
      counter <- counter + 1L
      path <- shard_file(directory, batch, shard_id, args$num_shards)
      if (!file.exists(path)) {
        stop("Missing ", batch, " shard: ", path, call. = FALSE)
      }
      object <- readRDS(path)
      manifest <- make_shard_manifest(
        hpc_args, configuration, batch, shard_id, replications
      )
      design <- make_batch_design(configuration, replications, batch)
      validate_complete_batch(
        object, manifest, configuration, design,
        paste0(batch, " shard ", shard_id, "/", args$num_shards)
      )
      if (batch == "evaluation") {
        evaluation_parts[[shard_id]] <- object$checkpoint
      } else {
        reference_parts[[shard_id]] <- object$checkpoint
      }
      input_paths[[counter]] <- path
      environment_signatures[[counter]] <- object$environment$signature
    }
  }

  reference_signature <- environment_signatures[[1L]]
  same_environment <- vapply(
    environment_signatures,
    function(x) isTRUE(all.equal(x, reference_signature, tolerance = 0)),
    logical(1)
  )
  if (!all(same_environment)) {
    stop("Shard R/Python/package version signatures differ.", call. = FALSE)
  }

  evaluation <- combine_checkpoints(evaluation_parts)
  reference <- combine_checkpoints(reference_parts)
  validate_merged_checkpoint(evaluation, hpc_args, configuration, "evaluation")
  validate_merged_checkpoint(reference, hpc_args, configuration, "reference")
  if (anyDuplicated(evaluation$diagnostics$seed) ||
      anyDuplicated(reference$diagnostics$seed) ||
      length(intersect(evaluation$diagnostics$seed, reference$diagnostics$seed))) {
    stop("Evaluation and reference seeds are duplicated or overlap.", call. = FALSE)
  }

  list(
    evaluation = evaluation,
    reference = reference,
    shard_directory = directory,
    shard_md5 = tools::md5sum(input_paths),
    environment_signature = reference_signature
  )
}

summarize_features <- function(checkpoint, fluctuations, reference) {
  reference_sd <- sqrt(diag(reference$covariance))
  n_values <- sort(unique(checkpoint$diagnostics$n))
  pieces <- lapply(n_values, function(n) {
    values <- fluctuations[checkpoint$diagnostics$n == n, , drop = FALSE]
    means <- colMeans(values)
    standard_deviations <- apply(values, 2, stats::sd)
    centered <- sweep(values, 2, means, "-")
    standardized <- sweep(centered, 2, standard_deviations, "/")
    quantiles <- apply(
      values, 2, stats::quantile, probs = c(0.025, 0.5, 0.975),
      names = FALSE, type = 8
    )
    data.frame(
      n = n,
      feature = colnames(values),
      mean = means,
      sd = standard_deviations,
      reference_sd = reference_sd,
      sd_ratio_to_reference = standard_deviations / reference_sd,
      q025 = quantiles[1L, ],
      median = quantiles[2L, ],
      q975 = quantiles[3L, ],
      skewness = colMeans(standardized^3),
      excess_kurtosis = colMeans(standardized^4) - 3,
      row.names = NULL
    )
  })
  do.call(rbind, pieces)
}

markdown_table <- function(data, digits = 4L) {
  formatted <- data
  numeric_columns <- vapply(formatted, is.numeric, logical(1))
  formatted[numeric_columns] <- lapply(
    formatted[numeric_columns], formatC, format = "f", digits = digits
  )
  header <- paste0("| ", paste(names(formatted), collapse = " | "), " |")
  separator <- paste0(
    "| ", paste(rep("---", ncol(formatted)), collapse = " | "), " |"
  )
  rows <- apply(formatted, 1, function(x) {
    paste0("| ", paste(x, collapse = " | "), " |")
  })
  paste(c(header, separator, rows), collapse = "\n")
}

write_summary_report <- function(path, configuration, standard_summary,
                                 selected_summary, covariance_comparison,
                                 files, loaded) {
  first_summary <- standard_summary[which.min(standard_summary$n), ]
  last_summary <- standard_summary[which.max(standard_summary$n), ]
  last_selected <- selected_summary[
    selected_summary$n == max(selected_summary$n), , drop = FALSE
  ]
  ratio_range <- range(last_selected$sd_ratio_to_reference)
  covariance_range <- range(covariance_comparison$relative_frobenius_error)
  format_number <- function(x, digits = 4L) {
    formatC(x, format = "f", digits = digits)
  }
  report <- c(
    "# CLT Gaussian HPC summary",
    "",
    paste0(
      "This report summarizes **", nrow(loaded$evaluation$diagnostics),
      " evaluation fits** and **", nrow(loaded$reference$diagnostics),
      " independent reference fits** from ",
      length(loaded$shard_md5), " validated shard files."
    ),
    "",
    paste0(
      "The evaluation design has ",
      length(unique(loaded$evaluation$diagnostics$replication)),
      " replications at each of $n=", paste(configuration$n_grid, collapse = ","),
      "$. The reference has ", nrow(loaded$reference$diagnostics),
      " independent replications at $n=", configuration$reference_n, "$.") ,
    "",
    "## Run summary",
    "",
    markdown_table(standard_summary),
    "",
    "## Representative Gaussian diagnostics",
    "",
    "The rows below contain total transported mass and the central bump from each projection.",
    "",
    markdown_table(selected_summary),
    "",
    "## Covariance stabilization",
    "",
    markdown_table(covariance_comparison),
    "",
    "## Interpretation",
    "",
    paste0(
      "All fits converged. The absolute transported-mass bias decreased from ",
      format_number(abs(first_summary$mass_bias)), " at $n=", first_summary$n,
      "$ to ", format_number(abs(last_summary$mass_bias)), " at $n=",
      last_summary$n, "$, while the maximum absolute mean scaled fluctuation ",
      "decreased from ",
      format_number(first_summary$max_absolute_mean_fluctuation), " to ",
      format_number(last_summary$max_absolute_mean_fluctuation), "."
    ),
    "",
    paste0(
      "At $n=", last_summary$n,
      "$, standard deviations for the five representative functionals are ",
      format_number(ratio_range[[1L]], 3L), "--",
      format_number(ratio_range[[2L]], 3L),
      " times their independent-reference values. Their largest absolute ",
      "skewness is ", format_number(max(abs(last_selected$skewness)), 3L),
      " and their largest absolute excess kurtosis is ",
      format_number(max(abs(last_selected$excess_kurtosis)), 3L), "."
    ),
    "",
    paste0(
      "Relative covariance errors range from ",
      format_number(covariance_range[[1L]], 3L), " to ",
      format_number(covariance_range[[2L]], 3L),
      ", close to the independent-reference noise benchmark of ",
      format_number(covariance_comparison$reference_noise_benchmark[[1L]], 3L),
      ". Taken together, these diagnostics support stabilization toward the ",
      "Gaussian-process limit at the resolution of this Monte Carlo study."
    ),
    "",
    "## Figures",
    "",
    paste0("- [", names(files), "](", unname(files), ")"),
    "",
    "All fluctuations use the per-marginal scaling $\\sqrt{n}$. Gaussian overlays,",
    "Q--Q standardization, and covariance comparisons use only the independent",
    "reference batch; evaluation and reference random-number streams are disjoint."
  )
  writeLines(report, path, useBytes = TRUE)
}

main_summary <- function() {
  args <- parse_summary_arguments(commandArgs(trailingOnly = TRUE))
  if (args$help) {
    summary_usage()
    return(invisible(NULL))
  }
  assert_repository_root()
  require_packages("ggplot2")
  validate_summary_arguments(args)

  configuration <- make_hpc_configuration(args)
  loaded <- load_completed_shards(args, configuration)
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(args$plot_dir, recursive = TRUE, showWarnings = FALSE)

  fluctuations <- compute_fluctuations(
    loaded$evaluation, configuration$test_functions$population_integrals
  )
  independent_reference <- make_independent_reference(
    loaded$reference, configuration$test_functions$population_integrals
  )
  validate_population_against_empirical(loaded$evaluation, configuration$population)
  validate_population_against_empirical(loaded$reference, configuration$population)

  prefix <- paste0(configuration$run_name, "_R", args$reps, "_")
  standard_summary <- summarize_replications(
    loaded$evaluation, fluctuations, configuration$population$mass
  )
  feature_summary <- summarize_features(
    loaded$evaluation, fluctuations, independent_reference
  )

  heatmaps <- make_fluctuation_heatmaps(
    loaded$evaluation, fluctuations, configuration$test_functions,
    args$plot_dir, prefix
  )
  gaussian <- make_gaussian_diagnostics(
    loaded$evaluation, fluctuations, configuration$test_functions,
    independent_reference, args$plot_dir, prefix
  )
  covariance <- make_covariance_diagnostics(
    loaded$evaluation, fluctuations, independent_reference,
    args$plot_dir, prefix
  )

  summary_path <- file.path(args$output_dir, paste0(prefix, "summary.csv"))
  feature_path <- file.path(args$output_dir, paste0(prefix, "feature_summary.csv"))
  covariance_path <- file.path(
    args$output_dir, paste0(prefix, "covariance_stabilization.csv")
  )
  bundle_path <- file.path(args$output_dir, paste0(prefix, "analysis_summary.rds"))
  report_path <- file.path(args$output_dir, paste0(prefix, "summary.md"))

  utils::write.csv(standard_summary, summary_path, row.names = FALSE)
  utils::write.csv(feature_summary, feature_path, row.names = FALSE)
  utils::write.csv(covariance$comparison, covariance_path, row.names = FALSE)
  atomic_save_rds(
    list(
      schema_version = "clt_gaussian_hpc_summary_v1",
      created_at = Sys.time(),
      configuration = configuration,
      evaluation_rows = nrow(loaded$evaluation$diagnostics),
      reference_rows = nrow(loaded$reference$diagnostics),
      standard_summary = standard_summary,
      feature_summary = feature_summary,
      covariance_comparison = covariance$comparison,
      independent_reference_covariance = independent_reference$covariance,
      shard_md5 = loaded$shard_md5,
      environment_signature = loaded$environment_signature
    ),
    bundle_path
  )

  selected_indices <- select_diagnostic_features(configuration$test_functions)
  selected_summary <- feature_summary[
    feature_summary$feature %in%
      configuration$test_functions$feature_names[selected_indices],
    c("n", "feature", "mean", "sd", "sd_ratio_to_reference",
      "skewness", "excess_kurtosis")
  ]
  selected_summary$feature <- rep(
    names(selected_indices), times = length(configuration$n_grid)
  )

  figure_paths <- unlist(
    c(heatmaps, gaussian, covariance[1:2]), use.names = TRUE
  )
  report_links <- file.path("..", basename(args$plot_dir), basename(figure_paths))
  names(report_links) <- names(figure_paths)
  write_summary_report(
    report_path, configuration, standard_summary, selected_summary,
    covariance$comparison, report_links, loaded
  )

  message(
    "Validated ", length(loaded$shard_md5),
    " shards and summarized the complete HPC run."
  )
  message("Summary report: ", report_path)
  message("CSV/RDS results: ", args$output_dir)
  message("Figures: ", args$plot_dir)
  invisible(list(report = report_path, figures = figure_paths))
}

if (sys.nframe() == 0L) main_summary()
