#!/usr/bin/env Rscript

# Combine the n = 10, 20 extension with the completed n = 50:500 run.

options(stringsAsFactors = FALSE)

summary_script <- file.path(
  "simulation", "ET_round2", "CLT_gaussian_summary.R"
)
small_n_worker_script <- file.path(
  "simulation", "ET_round2", "CLT_gaussian_small_n_hpc.R"
)
if (!all(file.exists(c(summary_script, small_n_worker_script)))) {
  stop("Run this script from the repository root.", call. = FALSE)
}
source(summary_script)
source(small_n_worker_script)

parse_combined_arguments <- function(args) {
  value_after_equals <- function(prefix, default) {
    hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
    if (!length(hit)) return(default)
    sub(paste0("^", prefix, "="), "", hit[[length(hit)]])
  }
  target_reps <- as.integer(value_after_equals("--reps", "10000"))
  base_reps <- as.integer(value_after_equals("--base-reps", "10000"))
  list(
    help = "--help" %in% args || "-h" %in% args,
    smoke_test = FALSE,
    reps = target_reps,
    base_reps = base_reps,
    base_num_shards = as.integer(
      value_after_equals("--base-num-shards", "100")
    ),
    extension_rep_start = as.integer(value_after_equals(
      "--extension-rep-start", as.character(base_reps + 1L)
    )),
    extension_num_shards = as.integer(
      value_after_equals("--extension-num-shards", "100")
    ),
    seed = as.integer(value_after_equals("--seed", "20260814")),
    shard_dir = value_after_equals(
      "--base-shard-dir",
      "/cwork/yx306/UOT/out/clt-gaussian-10000/shards"
    ),
    small_n_shard_dir = value_after_equals(
      "--small-n-shard-dir",
      "/cwork/yx306/UOT/out/clt-gaussian-10000/small-n-shards"
    ),
    output_dir = value_after_equals(
      "--output-dir", "/cwork/yx306/UOT/out/clt-gaussian-10000/result"
    ),
    plot_dir = value_after_equals(
      "--plot-dir", "/cwork/yx306/UOT/out/clt-gaussian-10000/plot"
    ),
    mc_block_size = as.integer(value_after_equals(
      "--mc-block-size", if (target_reps > 10000L) "2000" else "400"
    ))
  )
}

combined_usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript --vanilla simulation/ET_round2/CLT_gaussian_small_n_summary.R \\\n",
    "    --reps=100000 --base-reps=10000 \\\n",
    "    --base-num-shards=100 --extension-num-shards=100 \\\n",
    "    --extension-rep-start=10001 \\\n",
    "    --base-shard-dir=/cwork/yx306/UOT/out/clt-gaussian-10000/shards \\\n",
    "    --small-n-shard-dir=/cwork/yx306/UOT/out/clt-gaussian-10000/small-n-shards \\\n",
    "    --output-dir=/cwork/yx306/UOT/out/clt-gaussian-10000/result \\\n",
    "    --plot-dir=/cwork/yx306/UOT/out/clt-gaussian-10000/plot\n",
    sep = ""
  )
}

validate_combined_arguments <- function(args) {
  integers <- c(
    args$reps, args$base_reps, args$base_num_shards,
    args$extension_rep_start, args$extension_num_shards,
    args$seed, args$mc_block_size
  )
  if (anyNA(integers) || args$reps < 4L || args$base_reps < 4L ||
      args$base_reps > args$reps || args$base_num_shards < 1L ||
      args$extension_num_shards < 1L ||
      args$seed < 0L || args$mc_block_size < 4L ||
      args$reps %% args$mc_block_size != 0L) {
    stop(
      "Invalid arguments; --mc-block-size must divide --reps exactly.",
      call. = FALSE
    )
  }
  if (args$reps > args$base_reps &&
      args$extension_rep_start != args$base_reps + 1L) {
    stop("Continuation must start at --base-reps + 1.", call. = FALSE)
  }
  allowed_root <- "/cwork/yx306/UOT/out/"
  paths <- c(
    args$shard_dir, args$small_n_shard_dir,
    args$output_dir, args$plot_dir
  )
  normalized <- vapply(paths, normalizePath, character(1), mustWork = FALSE)
  if (any(!startsWith(normalized, allowed_root))) {
    stop("Every input and output path must be below /cwork/yx306/UOT/out.",
         call. = FALSE)
  }
  invisible(TRUE)
}

load_small_n_shards <- function(run_args, configuration,
                                expected_environment, label) {
  directory <- configuration_directory(run_args, configuration)
  if (!dir.exists(directory)) {
    stop("Small-n shard directory does not exist: ", directory, call. = FALSE)
  }

  parts <- vector("list", run_args$num_shards)
  paths <- character(run_args$num_shards)
  signatures <- vector("list", run_args$num_shards)
  for (shard_id in seq_len(run_args$num_shards)) {
    replications <- assigned_replications(run_args, shard_id)
    path <- shard_file(
      directory, "evaluation", shard_id, run_args$num_shards
    )
    if (!file.exists(path)) {
      stop("Missing small-n evaluation shard: ", path, call. = FALSE)
    }
    object <- readRDS(path)
    manifest <- make_small_n_manifest(
      run_args, configuration, shard_id, replications
    )
    design <- make_batch_design(configuration, replications, "evaluation")
    validate_complete_batch(
      object, manifest, configuration, design,
      paste0(label, " shard ", shard_id, "/", run_args$num_shards)
    )
    parts[[shard_id]] <- object$checkpoint
    paths[[shard_id]] <- path
    signatures[[shard_id]] <- object$environment$signature
  }

  same_environment <- vapply(
    signatures,
    function(x) isTRUE(all.equal(x, expected_environment, tolerance = 0)),
    logical(1)
  )
  if (!all(same_environment)) {
    stop("Small-n and original shards used different environments.",
         call. = FALSE)
  }

  checkpoint <- combine_checkpoints(parts)
  validate_checkpoint(
    checkpoint, configuration$test_functions$feature_names,
    configuration$n_grid, run_args$seed, 0L, label
  )
  expected <- make_batch_design(
    configuration, seq.int(run_args$rep_start, run_args$reps), "evaluation"
  )
  if (!setequal(job_keys(checkpoint$diagnostics), job_keys(expected)) ||
      nrow(checkpoint$diagnostics) != nrow(expected)) {
    stop("Small-n shards do not contain the exact requested design.",
         call. = FALSE)
  }
  list(
    checkpoint = checkpoint,
    paths = paths,
    md5 = tools::md5sum(paths),
    directory = directory
  )
}

load_hpc_extension_shards <- function(run_args, configuration,
                                      expected_environment) {
  directory <- configuration_directory(run_args, configuration)
  if (!dir.exists(directory)) {
    stop("Continuation shard directory does not exist: ", directory,
         call. = FALSE)
  }
  evaluation_parts <- vector("list", run_args$num_shards)
  reference_parts <- vector("list", run_args$num_shards)
  paths <- character(2L * run_args$num_shards)
  signatures <- vector("list", 2L * run_args$num_shards)
  counter <- 0L
  for (shard_id in seq_len(run_args$num_shards)) {
    replications <- assigned_replications(run_args, shard_id)
    for (batch in c("evaluation", "reference")) {
      counter <- counter + 1L
      path <- shard_file(directory, batch, shard_id, run_args$num_shards)
      if (!file.exists(path)) {
        stop("Missing continuation ", batch, " shard: ", path,
             call. = FALSE)
      }
      object <- readRDS(path)
      manifest <- make_shard_manifest(
        run_args, configuration, batch, shard_id, replications
      )
      design <- make_batch_design(configuration, replications, batch)
      validate_complete_batch(
        object, manifest, configuration, design,
        paste0(
          "continuation ", batch, " shard ", shard_id, "/",
          run_args$num_shards
        )
      )
      if (batch == "evaluation") {
        evaluation_parts[[shard_id]] <- object$checkpoint
      } else {
        reference_parts[[shard_id]] <- object$checkpoint
      }
      paths[[counter]] <- path
      signatures[[counter]] <- object$environment$signature
    }
  }
  same_environment <- vapply(
    signatures,
    function(x) isTRUE(all.equal(x, expected_environment, tolerance = 0)),
    logical(1)
  )
  if (!all(same_environment)) {
    stop("Continuation and base shards used different environments.",
         call. = FALSE)
  }
  replications <- seq.int(run_args$rep_start, run_args$reps)
  evaluation <- combine_checkpoints(evaluation_parts)
  reference <- combine_checkpoints(reference_parts)
  for (batch in c("evaluation", "reference")) {
    checkpoint <- if (batch == "evaluation") evaluation else reference
    allowed_n <- if (batch == "evaluation") {
      configuration$n_grid
    } else {
      configuration$reference_n
    }
    validate_checkpoint(
      checkpoint, configuration$test_functions$feature_names,
      allowed_n, run_args$seed, if (batch == "evaluation") 0L else 1L,
      paste0("Combined continuation ", batch)
    )
    expected <- make_batch_design(configuration, replications, batch)
    if (nrow(checkpoint$diagnostics) != nrow(expected) ||
        !setequal(job_keys(checkpoint$diagnostics), job_keys(expected))) {
      stop("Continuation ", batch, " design is incomplete.", call. = FALSE)
    }
  }
  list(
    evaluation = evaluation,
    reference = reference,
    paths = paths,
    md5 = tools::md5sum(paths),
    directory = directory
  )
}

calibrate_mc_uncertainty <- function(checkpoint, fluctuations, reference,
                                     block_size) {
  reference_n <- reference$n
  evaluation_values <- fluctuations[
    checkpoint$diagnostics$n == reference_n, , drop = FALSE
  ]
  reference_values <- reference$fluctuations
  target_replications <- nrow(reference_values)
  if (nrow(evaluation_values) != target_replications ||
      target_replications %% block_size != 0L) {
    stop("The n=500 null batches are incompatible with the block calibration.",
         call. = FALSE)
  }

  null_values <- rbind(evaluation_values, reference_values)
  block <- rep(
    seq_len(nrow(null_values) / block_size), each = block_size
  )
  block_covariances <- lapply(
    split(seq_len(nrow(null_values)), block),
    function(i) stats::cov(null_values[i, , drop = FALSE])
  )
  pairs <- utils::combn(seq_along(block_covariances), 2L)
  scale_to_target <- sqrt(
    (block_size - 1) / (target_replications - 1)
  )
  reference_norm <- norm(reference$covariance, type = "F")
  distances <- apply(pairs, 2L, function(pair) {
    norm(
      block_covariances[[pair[[1L]]]] -
        block_covariances[[pair[[2L]]]],
      type = "F"
    ) / reference_norm * scale_to_target
  })
  interval <- stats::quantile(
    distances, probs = c(0.025, 0.5, 0.975), names = FALSE, type = 8
  )
  rms <- sqrt(
    2 * (
      sum(reference$covariance^2) +
        sum(diag(reference$covariance))^2
    ) / (
      (target_replications - 1) * sum(reference$covariance^2)
    )
  )
  list(
    lower = interval[[1L]],
    median = interval[[2L]],
    upper = interval[[3L]],
    rms = rms,
    block_size = block_size,
    blocks = length(block_covariances),
    pairwise_distances = distances,
    target_replications = target_replications,
    method = paste0(
      "All pairwise covariance distances among ",
      length(block_covariances), " independent n=500 blocks of size ",
      block_size, ", scaled by sqrt((", block_size, "-1)/(",
      target_replications, "-1))"
    )
  )
}

make_frobenius_plot <- function(checkpoint, fluctuations, reference,
                                calibration, plot_dir, prefix) {
  n_values <- sort(unique(checkpoint$diagnostics$n))
  reference_norm <- norm(reference$covariance, type = "F")
  comparison <- do.call(rbind, lapply(n_values, function(n) {
    indices <- checkpoint$diagnostics$n == n
    covariance <- stats::cov(fluctuations[indices, , drop = FALSE])
    distance <- norm(
      covariance - reference$covariance, type = "F"
    ) / reference_norm
    data.frame(
      n = n,
      replications = sum(indices),
      relative_frobenius_distance = distance,
      mc_rms_benchmark = calibration$rms,
      mc_band_lower = calibration$lower,
      mc_band_upper = calibration$upper,
      mc_classification = if (distance > calibration$upper) {
        "Above MC band"
      } else {
        "Within MC band"
      }
    )
  }))
  comparison$mc_classification <- factor(
    comparison$mc_classification,
    levels = c("Above MC band", "Within MC band")
  )

  band <- data.frame(
    n = range(n_values),
    lower = calibration$lower,
    upper = calibration$upper
  )
  plot <- ggplot2::ggplot(
    comparison,
    ggplot2::aes(x = n, y = relative_frobenius_distance)
  ) +
    ggplot2::geom_ribbon(
      data = band,
      ggplot2::aes(x = n, ymin = lower, ymax = upper),
      inherit.aes = FALSE, fill = "#B8C7D9", alpha = 0.55
    ) +
    ggplot2::geom_hline(
      yintercept = calibration$rms,
      linetype = "dashed", color = "#8B1E2D", linewidth = 0.8
    ) +
    ggplot2::geom_line(color = "#30343B", linewidth = 0.8) +
    ggplot2::geom_point(
      ggplot2::aes(color = mc_classification), size = 2.8
    ) +
    ggplot2::scale_color_manual(
      values = c("Above MC band" = "#B2182B", "Within MC band" = "#2166AC"),
      drop = FALSE
    ) +
    ggplot2::scale_x_log10(breaks = n_values, labels = n_values) +
    ggplot2::labs(
      title = "Finite-sample Covariance Distance from the n = 500 Reference",
      subtitle = paste0(
        "Shaded: 95% Monte Carlo null band; dashed: RMS benchmark (R = ",
        calibration$target_replications, ")"
      ),
      x = "Sample size per marginal (n, log scale)",
      y = "Relative Frobenius distance",
      color = NULL,
      caption = paste0(
        "Null band: ", calibration$blocks,
        " independent n=500 blocks of size ", calibration$block_size,
        ", scaled to R = ", calibration$target_replications,
        ". The n=500 evaluation point is an independent noise-only control."
      )
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank()
    )

  path <- file.path(plot_dir, paste0(prefix, "frobenius_distance_mc_band.pdf"))
  ggplot2::ggsave(path, plot, width = 8.2, height = 5.6, units = "in")
  list(path = path, comparison = comparison)
}

write_combined_report <- function(path, summary, comparison, calibration,
                                  figure_path, replications) {
  above <- comparison$n[
    comparison$relative_frobenius_distance > calibration$upper
  ]
  report <- c(
    "# Small-n Gaussian UOT covariance convergence",
    "",
    paste0(
      "This analysis combines ", format(replications, big.mark = ","),
      " evaluation replications at each of"
    ),
    paste0(
      "$n\\in\\{10,20,50,100,200,500\\}$ with an independent ",
      format(replications, big.mark = ","), "-run reference batch at $n=500$."
    ),
    "",
    "## Run summary",
    "",
    markdown_table(summary),
    "",
    "## Relative Frobenius distance",
    "",
    markdown_table(comparison[, c(
      "n", "relative_frobenius_distance", "mc_rms_benchmark",
      "mc_band_lower", "mc_band_upper", "mc_classification"
    )]),
    "",
    paste0(
      "The 95% Monte Carlo null band is [",
      formatC(calibration$lower, format = "f", digits = 4), ", ",
      formatC(calibration$upper, format = "f", digits = 4),
      "], with RMS benchmark ",
      formatC(calibration$rms, format = "f", digits = 4), "."
    ),
    "",
    if (length(above)) {
      paste0(
        "Sample sizes above the calibrated Monte Carlo band: ",
        paste(above, collapse = ", "), "."
      )
    } else {
      "No sample size is above the calibrated Monte Carlo band."
    },
    "",
    paste0("[Frobenius-distance figure](../plot/", basename(figure_path), ")"),
    "",
    paste0(
      "Band construction: ", calibration$method,
      ". Pairwise block distances are dependent, so the band is an empirical ",
      "diagnostic rather than a formal confidence interval."
    )
  )
  writeLines(report, path, useBytes = TRUE)
}

main_combined_summary <- function() {
  args <- parse_combined_arguments(commandArgs(trailingOnly = TRUE))
  if (args$help) {
    combined_usage()
    return(invisible(NULL))
  }
  assert_repository_root()
  require_packages("ggplot2")
  validate_combined_arguments(args)

  base_args <- list(
    smoke_test = FALSE,
    reps = args$base_reps,
    rep_start = 1L,
    num_shards = args$base_num_shards,
    seed = args$seed,
    shard_dir = args$shard_dir
  )
  base_configuration <- make_hpc_configuration(base_args)
  base <- load_completed_shards(base_args, base_configuration)
  small_configuration <- make_small_n_configuration()
  base_small_args <- list(
    smoke_test = FALSE,
    reps = args$base_reps,
    rep_start = 1L,
    num_shards = args$base_num_shards,
    seed = args$seed,
    shard_dir = args$small_n_shard_dir
  )
  base_small <- load_small_n_shards(
    base_small_args, small_configuration, base$environment_signature,
    "base small-n evaluation"
  )
  evaluation_parts <- list(base_small$checkpoint, base$evaluation)
  reference_parts <- list(base$reference)
  continuation_md5 <- character()
  continuation_small_md5 <- character()
  if (args$reps > args$base_reps) {
    extension_args <- list(
      smoke_test = FALSE,
      reps = args$reps,
      rep_start = args$extension_rep_start,
      num_shards = args$extension_num_shards,
      seed = args$seed,
      shard_dir = args$shard_dir
    )
    continuation <- load_hpc_extension_shards(
      extension_args, base_configuration, base$environment_signature
    )
    extension_small_args <- extension_args
    extension_small_args$shard_dir <- args$small_n_shard_dir
    continuation_small <- load_small_n_shards(
      extension_small_args, small_configuration, base$environment_signature,
      "continuation small-n evaluation"
    )
    evaluation_parts <- c(
      evaluation_parts,
      list(continuation_small$checkpoint, continuation$evaluation)
    )
    reference_parts <- c(reference_parts, list(continuation$reference))
    continuation_md5 <- continuation$md5
    continuation_small_md5 <- continuation_small$md5
  }
  evaluation <- combine_checkpoints(evaluation_parts)
  reference_checkpoint <- combine_checkpoints(reference_parts)
  n_grid <- sort(unique(c(
    small_configuration$n_grid, base_configuration$n_grid
  )))
  validate_checkpoint(
    evaluation, base_configuration$test_functions$feature_names,
    n_grid, args$seed, 0L, "Combined n=10:500 evaluation"
  )
  expected <- expand.grid(
    n = n_grid, replication = seq_len(args$reps),
    KEEP.OUT.ATTRS = FALSE
  )
  if (nrow(evaluation$diagnostics) != nrow(expected) ||
      !setequal(job_keys(evaluation$diagnostics), job_keys(expected))) {
    stop("Combined evaluation does not contain the exact n=10:500 design.",
         call. = FALSE)
  }
  if (nrow(reference_checkpoint$diagnostics) != args$reps ||
      !setequal(
        reference_checkpoint$diagnostics$replication, seq_len(args$reps)
      ) ||
      any(reference_checkpoint$diagnostics$n != 500L)) {
    stop("Combined reference does not contain the exact requested design.",
         call. = FALSE)
  }
  if (anyDuplicated(evaluation$diagnostics$seed) ||
      anyDuplicated(reference_checkpoint$diagnostics$seed) ||
      length(intersect(
        evaluation$diagnostics$seed, reference_checkpoint$diagnostics$seed
      ))) {
    stop("Combined evaluation/reference seeds duplicate or overlap.",
         call. = FALSE)
  }

  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(args$plot_dir, recursive = TRUE, showWarnings = FALSE)
  fluctuations <- compute_fluctuations(
    evaluation, base_configuration$test_functions$population_integrals
  )
  reference <- make_independent_reference(
    reference_checkpoint,
    base_configuration$test_functions$population_integrals
  )
  calibration <- calibrate_mc_uncertainty(
    evaluation, fluctuations, reference, args$mc_block_size
  )
  prefix <- paste0("clt_gaussian_n10_to_500_R", args$reps, "_")
  figure <- make_frobenius_plot(
    evaluation, fluctuations, reference, calibration,
    args$plot_dir, prefix
  )
  summary <- summarize_replications(
    evaluation, fluctuations, base_configuration$population$mass
  )

  summary_path <- file.path(args$output_dir, paste0(prefix, "summary.csv"))
  comparison_path <- file.path(
    args$output_dir, paste0(prefix, "frobenius_distance_mc_band.csv")
  )
  calibration_path <- file.path(
    args$output_dir, paste0(prefix, "mc_calibration.rds")
  )
  report_path <- file.path(args$output_dir, paste0(prefix, "summary.md"))
  utils::write.csv(summary, summary_path, row.names = FALSE)
  utils::write.csv(figure$comparison, comparison_path, row.names = FALSE)
  atomic_save_rds(
    list(
      schema_version = "clt_gaussian_small_n_mc_calibration_v1",
      created_at = Sys.time(),
      calibration = calibration,
      base_shard_md5 = base$shard_md5,
      base_small_n_shard_md5 = base_small$md5,
      continuation_shard_md5 = continuation_md5,
      continuation_small_n_shard_md5 = continuation_small_md5,
      environment_signature = base$environment_signature
    ),
    calibration_path
  )
  write_combined_report(
    report_path, summary, figure$comparison, calibration, figure$path,
    args$reps
  )
  message(
    "Validated and combined ", format(args$reps, big.mark = ","),
    " replications at n = 10:500."
  )
  message("Report: ", report_path)
  message("Frobenius-distance figure: ", figure$path)
  invisible(list(report = report_path, figure = figure$path))
}

if (sys.nframe() == 0L) main_combined_summary()
