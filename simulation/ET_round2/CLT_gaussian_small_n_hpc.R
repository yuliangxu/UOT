#!/usr/bin/env Rscript

# Evaluation-only HPC extension for n = 10 and n = 20.

options(stringsAsFactors = FALSE)

core_script <- file.path("simulation", "ET_round2", "CLT_gaussian_experiment.R")
hpc_script <- file.path("simulation", "ET_round2", "CLT_gaussian_hpc.R")
small_n_script <- file.path(
  "simulation", "ET_round2", "CLT_gaussian_small_n_hpc.R"
)
if (!all(file.exists(c(core_script, hpc_script, small_n_script)))) {
  stop("Run this script from the repository root.", call. = FALSE)
}
source(core_script)
source(hpc_script)

parse_small_n_arguments <- function(args) {
  value_after_equals <- function(prefix, default) {
    hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
    if (!length(hit)) return(default)
    sub(paste0("^", prefix, "="), "", hit[[length(hit)]])
  }
  list(
    help = "--help" %in% args || "-h" %in% args,
    smoke_test = "--smoke-test" %in% args,
    mode = value_after_equals("--mode", ""),
    reps = as.integer(value_after_equals("--reps", "10000")),
    rep_start = as.integer(value_after_equals("--rep-start", "1")),
    num_shards = as.integer(value_after_equals("--num-shards", "100")),
    shard_id = as.integer(value_after_equals("--shard-id", NA_character_)),
    seed = as.integer(value_after_equals("--seed", "20260814")),
    shard_dir = value_after_equals(
      "--shard-dir",
      "/cwork/yx306/UOT/out/clt-gaussian-10000/small-n-shards"
    )
  )
}

small_n_usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript --vanilla simulation/ET_round2/CLT_gaussian_small_n_hpc.R \\\n",
    "    --mode=preflight\n",
    "  Rscript --vanilla simulation/ET_round2/CLT_gaussian_small_n_hpc.R \\\n",
    "    --mode=shard --reps=10000 --rep-start=1 \\\n",
    "    --num-shards=100 --shard-id=1 \\\n",
    "    --shard-dir=/cwork/yx306/UOT/out/clt-gaussian-10000/small-n-shards\n",
    sep = ""
  )
}

validate_small_n_arguments <- function(args) {
  if (!args$mode %in% c("preflight", "shard")) {
    stop("--mode must be preflight or shard.", call. = FALSE)
  }
  integers <- c(args$reps, args$rep_start, args$num_shards, args$seed)
  if (anyNA(integers) || args$reps < 1L || args$rep_start < 1L ||
      args$rep_start > args$reps || args$num_shards < 1L || args$seed < 0L) {
    stop("Invalid replication, shard, or seed arguments.", call. = FALSE)
  }
  number_to_assign <- args$reps - args$rep_start + 1L
  if (args$num_shards > number_to_assign) {
    stop("--num-shards exceeds the number of replications.", call. = FALSE)
  }
  if (args$mode == "shard" &&
      (is.na(args$shard_id) || args$shard_id < 1L ||
       args$shard_id > args$num_shards)) {
    stop("Shard mode requires 1 <= --shard-id <= --num-shards.", call. = FALSE)
  }
  normalized <- normalizePath(args$shard_dir, mustWork = FALSE)
  allowed <- startsWith(normalized, "/cwork/yx306/UOT/out/")
  if (!allowed && !args$smoke_test) {
    stop("Production shards must be below /cwork/yx306/UOT/out.", call. = FALSE)
  }
  invisible(TRUE)
}

make_small_n_configuration <- function() {
  marginals <- make_marginals()
  population <- gaussian_uot_population(
    marginals$a, marginals$b, marginals$A, marginals$B, 0.5, 1.0
  )
  list(
    epsilon = 0.5,
    rho = 1.0,
    n_grid = c(10L, 20L),
    reference_n = 500L,
    grid_size = 7L,
    run_name = "clt_gaussian_small_n",
    marginals = marginals,
    population = population,
    test_functions = make_test_functions(population, 7L)
  )
}

small_n_source_provenance <- function() {
  provenance <- source_provenance()
  provenance$extension_source_md5 <- as.character(
    tools::md5sum(small_n_script)
  )
  provenance
}

make_small_n_manifest <- function(args, configuration, shard_id,
                                  replications) {
  manifest <- make_shard_manifest(
    args, configuration, "evaluation", shard_id, replications
  )
  manifest$schema_version <- "clt_small_n_hpc_shard_v1"
  manifest$source_provenance <- small_n_source_provenance()
  manifest
}

run_small_n_shard <- function(args, configuration) {
  require_packages(c("MASS", "ggplot2", "reticulate"))
  directory <- configuration_directory(args, configuration)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  replications <- assigned_replications(args, args$shard_id)
  environment <- environment_record()
  validate_worker_environment(environment)

  seeds <- unlist(lapply(configuration$n_grid, function(n) {
    vapply(
      replications,
      function(replication) cell_seed(args$seed, n, replication, 0L),
      integer(1)
    )
  }))
  if (anyDuplicated(seeds)) {
    stop("The planned small-n shard contains duplicate seeds.", call. = FALSE)
  }

  path <- shard_file(
    directory, "evaluation", args$shard_id, args$num_shards
  )
  manifest <- make_small_n_manifest(
    args, configuration, args$shard_id, replications
  )
  design <- make_batch_design(configuration, replications, "evaluation")
  label <- paste0(
    "small-n evaluation shard ", args$shard_id, "/", args$num_shards
  )
  if (file.exists(path)) {
    validate_complete_batch(
      readRDS(path), manifest, configuration, design, label
    )
    message("Validated existing ", label, ": ", path)
    return(invisible(path))
  }

  checkpoint <- empty_checkpoint(configuration$test_functions$feature_names)
  message("Running ", label, " with ", nrow(design), " jobs.")
  for (row in seq_len(nrow(design))) {
    job <- design[row, ]
    if (row == 1L || row %% 25L == 0L || row == nrow(design)) {
      message(
        "[", label, " ", row, "/", nrow(design), "] n=", job$n,
        ", replication=", job$replication
      )
    }
    result <- run_one_replication(
      job$n, job$replication, args$seed,
      configuration$marginals, configuration$population,
      configuration$test_functions, configuration$epsilon,
      configuration$rho, 0L
    )
    checkpoint <- append_result(checkpoint, result)
  }
  object <- list(
    manifest = manifest,
    environment = environment,
    checkpoint = checkpoint
  )
  validate_complete_batch(object, manifest, configuration, design, label)
  atomic_save_rds(object, path)
  message("Saved ", label, ": ", path)
  invisible(path)
}

run_small_n_preflight <- function(args, configuration) {
  require_packages(c("MASS", "ggplot2", "reticulate"))
  dir.create(args$shard_dir, recursive = TRUE, showWarnings = FALSE)
  environment <- environment_record()
  validate_worker_environment(environment)
  result <- run_one_replication(
    10L, 1L, args$seed, configuration$marginals,
    configuration$population, configuration$test_functions,
    configuration$epsilon, configuration$rho, 0L
  )
  if (!isTRUE(result$diagnostics$converged) ||
      result$diagnostics$final_error > 1e-9 ||
      any(!is.finite(result$integrals))) {
    stop("Small-n preflight solve failed.", call. = FALSE)
  }
  cat("Small-n HPC preflight passed for n = 10, 20.\n")
  cat("Python:", environment$Python_executable, "\n")
  invisible(environment)
}

main_small_n_hpc <- function() {
  args <- parse_small_n_arguments(commandArgs(trailingOnly = TRUE))
  if (args$help) {
    small_n_usage()
    return(invisible(NULL))
  }
  assert_repository_root()
  validate_small_n_arguments(args)
  configuration <- make_small_n_configuration()
  if (args$mode == "preflight") {
    run_small_n_preflight(args, configuration)
  } else {
    run_small_n_shard(args, configuration)
  }
  invisible(NULL)
}

if (sys.nframe() == 0L) main_small_n_hpc()
