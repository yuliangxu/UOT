#!/usr/bin/env Rscript

# Scheduler-neutral shard/merge runner for CLT_gaussian_experiment.R.
#
# Each array task writes files unique to its shard. No worker writes the
# canonical checkpoints, summaries, covariance, or plots. After every shard
# is complete, merge mode validates the exact global design, writes the
# canonical checkpoints, and invokes the ordinary experiment analysis to
# regenerate all final diagnostics.

options(stringsAsFactors = FALSE)

core_script <- file.path(
  "simulation", "ET_round2", "CLT_gaussian_experiment.R"
)
hpc_script <- file.path(
  "simulation", "ET_round2", "CLT_gaussian_hpc.R"
)
if (!file.exists(core_script)) {
  stop(
    "Run this script from the repository root. Expected to find ",
    core_script, call. = FALSE
  )
}
source(core_script)

source_provenance <- function() {
  paths <- c(core = core_script, hpc = hpc_script)
  hashes <- as.character(tools::md5sum(paths))
  names(hashes) <- names(paths)
  git_commit <- tryCatch(
    {
      value <- system2(
        "git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE
      )
      if (!is.null(attr(value, "status")) || !length(value)) {
        NA_character_
      } else {
        value[[1L]]
      }
    },
    error = function(e) NA_character_,
    warning = function(w) NA_character_
  )
  list(source_md5 = hashes, git_commit = git_commit)
}

hpc_usage <- function() {
  cat(
    paste0(
      "Usage:\n",
      "  Rscript simulation/ET_round2/CLT_gaussian_hpc.R --mode=preflight\n",
      "  Rscript simulation/ET_round2/CLT_gaussian_hpc.R --mode=shard \\\n",
      "    --reps=10000 --rep-start=1 --num-shards=100 --shard-id=1 \\\n",
      "    --shard-dir=/shared/path/clt_shards\n",
      "  Rscript simulation/ET_round2/CLT_gaussian_hpc.R --mode=merge \\\n",
      "    --reps=10000 --rep-start=1 --num-shards=100 \\\n",
      "    --shard-dir=/shared/path/clt_shards \\\n",
      "    --output-dir=/shared/path/clt_result\n\n",
      "Use --rep-start=2001 only when canonical 1:2000 evaluation and ",
      "reference checkpoints have been copied to --output-dir.\n"
    )
  )
}

parse_hpc_arguments <- function(args) {
  value_after_equals <- function(prefix, default = NULL) {
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
    num_shards = as.integer(value_after_equals("--num-shards", NA_character_)),
    shard_id = as.integer(value_after_equals("--shard-id", NA_character_)),
    seed = as.integer(value_after_equals("--seed", "20260814")),
    shard_dir = value_after_equals(
      "--shard-dir", "simulation/ET_round2/result/clt_shards"
    ),
    output_dir = value_after_equals(
      "--output-dir", "simulation/ET_round2/result"
    ),
    plot_dir = value_after_equals(
      "--plot-dir", "simulation/ET_round2/plot"
    )
  )
}

validate_hpc_arguments <- function(args) {
  if (args$help) return(invisible(TRUE))
  if (!args$mode %in% c("preflight", "shard", "merge")) {
    stop("--mode must be preflight, shard, or merge.", call. = FALSE)
  }
  if (is.na(args$seed) || args$seed < 0L ||
      args$seed >= .Machine$integer.max) {
    stop(
      "--seed must be an integer between 0 and 2147483646.",
      call. = FALSE
    )
  }
  if (args$mode == "preflight") return(invisible(TRUE))
  integer_arguments <- c(args$reps, args$rep_start, args$num_shards)
  if (anyNA(integer_arguments) || any(integer_arguments < 1L)) {
    stop("--reps, --rep-start, and --num-shards must be positive integers.",
         call. = FALSE)
  }
  if (args$rep_start > args$reps) {
    stop("--rep-start cannot exceed --reps.", call. = FALSE)
  }
  if (args$smoke_test && args$reps != 4L) {
    stop(
      "Smoke mode uses exactly four replications; set --reps=4.",
      call. = FALSE
    )
  }
  number_to_assign <- args$reps - args$rep_start + 1L
  if (args$num_shards > number_to_assign) {
    stop("--num-shards cannot exceed the number of assigned replications.",
         call. = FALSE)
  }
  if (args$mode == "shard" &&
      (is.na(args$shard_id) || args$shard_id < 1L ||
       args$shard_id > args$num_shards)) {
    stop("Shard mode requires 1 <= --shard-id <= --num-shards.",
         call. = FALSE)
  }
  invisible(TRUE)
}

make_hpc_configuration <- function(args) {
  epsilon <- 0.5
  rho <- 1.0
  if (args$smoke_test) {
    n_grid <- c(20L, 30L)
    grid_size <- 5L
    run_name <- "clt_gaussian_all_projections_smoke_test"
  } else {
    n_grid <- c(50L, 100L, 200L, 500L)
    grid_size <- 7L
    run_name <- "clt_gaussian_all_projections"
  }
  marginals <- make_marginals()
  population <- gaussian_uot_population(
    marginals$a, marginals$b, marginals$A, marginals$B, epsilon, rho
  )
  test_functions <- make_test_functions(population, grid_size)
  list(
    epsilon = epsilon,
    rho = rho,
    n_grid = n_grid,
    reference_n = max(n_grid),
    grid_size = grid_size,
    run_name = run_name,
    marginals = marginals,
    population = population,
    test_functions = test_functions
  )
}

assigned_replications <- function(args, shard_id) {
  ids <- seq.int(args$rep_start, args$reps)
  ids[((ids - args$rep_start) %% args$num_shards) == (shard_id - 1L)]
}

configuration_directory <- function(args, configuration) {
  file.path(
    args$shard_dir,
    sprintf(
      "%s_R%05d_from%05d_S%04d",
      configuration$run_name, args$reps, args$rep_start, args$num_shards
    )
  )
}

shard_file <- function(directory, batch, shard_id, num_shards) {
  file.path(
    directory,
    sprintf("%s_shard_%04d_of_%04d.rds", batch, shard_id, num_shards)
  )
}

make_batch_design <- function(configuration, replications, batch) {
  if (batch == "evaluation") {
    expand.grid(
      n = configuration$n_grid,
      replication = replications,
      KEEP.OUT.ATTRS = FALSE
    )
  } else if (batch == "reference") {
    data.frame(
      n = rep(configuration$reference_n, length(replications)),
      replication = replications
    )
  } else {
    stop("Unknown batch: ", batch, call. = FALSE)
  }
}

job_keys <- function(diagnostics) {
  paste(diagnostics$n, diagnostics$replication, sep = ":")
}

environment_record <- function() {
  package_names <- c("MASS", "ggplot2", "reticulate")
  package_versions <- vapply(
    package_names, function(x) as.character(utils::packageVersion(x)),
    character(1)
  )
  python_configuration <- reticulate::py_config()
  ot <- reticulate::import("ot", convert = TRUE)
  numpy <- reticulate::import("numpy", convert = TRUE)
  scipy <- reticulate::import("scipy", convert = TRUE)
  thread_variables <- c(
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
  )
  requested_python <- Sys.getenv("RETICULATE_PYTHON", unset = NA_character_)
  thread_settings <- stats::setNames(
    Sys.getenv(thread_variables, unset = NA_character_), thread_variables
  )
  session <- utils::sessionInfo()
  list(
    signature = list(
      R_version = R.version.string,
      platform = R.version$platform,
      package_versions = package_versions,
      BLAS = session$BLAS,
      LAPACK = session$LAPACK,
      Python_executable = python_configuration$python,
      Python_library = python_configuration$libpython,
      Python_version = as.character(python_configuration$version_string),
      POT_version = as.character(ot[["__version__"]]),
      NumPy_version = as.character(numpy[["__version__"]]),
      SciPy_version = as.character(scipy[["__version__"]]),
      RNG_kind = RNGkind(),
      requested_python = requested_python,
      thread_settings = thread_settings
    ),
    Python_executable = python_configuration$python,
    reticulate_requested_python = requested_python,
    thread_settings = thread_settings,
    recorded_at = Sys.time()
  )
}

validate_worker_environment <- function(record) {
  expected_rng <- c("Mersenne-Twister", "Inversion", "Rejection")
  if (!identical(record$signature$RNG_kind, expected_rng)) {
    stop(
      "Unexpected RNG configuration: ",
      paste(record$signature$RNG_kind, collapse = ", "),
      call. = FALSE
    )
  }

  requested <- record$reticulate_requested_python
  if (is.na(requested) || !nzchar(requested) ||
      !startsWith(path.expand(requested), "/")) {
    stop(
      "Set RETICULATE_PYTHON to an absolute Python executable path before ",
      "running preflight or a shard.", call. = FALSE
    )
  }
  requested <- normalizePath(path.expand(requested), mustWork = TRUE)
  selected <- normalizePath(record$Python_executable, mustWork = TRUE)
  if (!identical(requested, selected)) {
    stop(
      "reticulate selected ", selected,
      " instead of RETICULATE_PYTHON=", requested, ".",
      call. = FALSE
    )
  }

  values <- record$thread_settings
  if (anyNA(values) || any(values != "1")) {
    stop(
      "Set OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, MKL_NUM_THREADS, ",
      "VECLIB_MAXIMUM_THREADS, and NUMEXPR_NUM_THREADS to 1 before R starts.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

make_shard_manifest <- function(args, configuration, batch, shard_id,
                                replications) {
  list(
    schema_version = "clt_hpc_shard_v1",
    run_name = configuration$run_name,
    smoke_test = args$smoke_test,
    batch = batch,
    target_replications = args$reps,
    rep_start = args$rep_start,
    num_shards = args$num_shards,
    shard_id = shard_id,
    assigned_replications = as.integer(replications),
    master_seed = args$seed,
    seed_formula_version = "cell_seed_v2_stream_offset_1000000007",
    seed_stream = if (batch == "evaluation") 0L else 1L,
    n_values = if (batch == "evaluation") {
      configuration$n_grid
    } else {
      configuration$reference_n
    },
    epsilon = configuration$epsilon,
    rho = configuration$rho,
    cost = "fixed squared Euclidean distance",
    solver = list(
      method = "sinkhorn_stabilized", reg_type = "kl",
      max_iter = 10000L, tolerance = 1e-10
    ),
    fluctuation_scaling = "sqrt(n)",
    marginals = configuration$marginals,
    population_mass = configuration$population$mass,
    population_mean = configuration$population$mean,
    population_covariance = configuration$population$covariance,
    feature_names = configuration$test_functions$feature_names,
    population_integrals =
      configuration$test_functions$population_integrals,
    source_provenance = source_provenance()
  )
}

manifests_match <- function(observed, expected) {
  isTRUE(all.equal(observed, expected, tolerance = 1e-12))
}

validate_complete_batch <- function(object, expected_manifest,
                                    configuration, expected_design, label) {
  if (!is.list(object) || is.null(object$manifest) ||
      is.null(object$checkpoint) || is.null(object$environment)) {
    stop(label, " has an invalid shard structure.", call. = FALSE)
  }
  if (!manifests_match(object$manifest, expected_manifest)) {
    stop(label, " manifest does not match the requested configuration.",
         call. = FALSE)
  }
  checkpoint <- object$checkpoint
  validate_checkpoint(
    checkpoint,
    configuration$test_functions$feature_names,
    expected_manifest$n_values,
    expected_manifest$master_seed,
    expected_manifest$seed_stream,
    label
  )
  expected_keys <- job_keys(expected_design)
  observed_keys <- job_keys(checkpoint$diagnostics)
  if (length(observed_keys) != length(expected_keys) ||
      !setequal(observed_keys, expected_keys)) {
    stop(label, " does not contain its exact assigned job set.", call. = FALSE)
  }
  if (!all(is.finite(checkpoint$integrals)) ||
      !all(is.finite(checkpoint$diagnostics$plan_mass)) ||
      any(checkpoint$diagnostics$plan_mass <= 0)) {
    stop(label, " contains non-finite integrals or invalid plan masses.",
         call. = FALSE)
  }
  if (!isTRUE(all.equal(
    as.numeric(checkpoint$integrals[, "total_mass"]),
    checkpoint$diagnostics$plan_mass,
    tolerance = 0
  ))) {
    stop(label, " total-mass integrals do not align with diagnostics.",
         call. = FALSE)
  }
  if (!all(checkpoint$diagnostics$converged) ||
      any(!is.finite(checkpoint$diagnostics$final_error)) ||
      any(checkpoint$diagnostics$final_error > 1e-9)) {
    stop(label, " contains a nonconverged solver result.", call. = FALSE)
  }
  invisible(TRUE)
}

run_batch_shard <- function(args, configuration, environment, batch,
                            shard_id, replications, directory) {
  path <- shard_file(directory, batch, shard_id, args$num_shards)
  manifest <- make_shard_manifest(
    args, configuration, batch, shard_id, replications
  )
  design <- make_batch_design(configuration, replications, batch)
  label <- paste0(batch, " shard ", shard_id, "/", args$num_shards)
  if (file.exists(path)) {
    existing <- readRDS(path)
    validate_complete_batch(
      existing, manifest, configuration, design, label
    )
    message("Validated existing ", label, ": ", path)
    return(invisible(path))
  }

  checkpoint <- empty_checkpoint(
    configuration$test_functions$feature_names
  )
  seed_stream <- if (batch == "evaluation") 0L else 1L
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
      configuration$rho, seed_stream
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

run_shard_mode <- function(args, configuration) {
  require_packages(c("MASS", "ggplot2", "reticulate"))
  directory <- configuration_directory(args, configuration)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  replications <- assigned_replications(args, args$shard_id)
  environment <- environment_record()
  validate_worker_environment(environment)

  evaluation_seeds <- vapply(
    configuration$n_grid,
    function(n) {
      vapply(
        replications,
        function(replication) cell_seed(args$seed, n, replication, 0L),
        integer(1)
      )
    },
    integer(length(replications))
  )
  reference_seeds <- vapply(
    replications,
    function(replication) {
      cell_seed(args$seed, configuration$reference_n, replication, 1L)
    },
    integer(1)
  )
  if (anyDuplicated(c(evaluation_seeds)) || anyDuplicated(reference_seeds) ||
      length(intersect(c(evaluation_seeds), reference_seeds))) {
    stop("The planned shard contains duplicate or overlapping seeds.",
         call. = FALSE)
  }

  run_batch_shard(
    args, configuration, environment, "evaluation", args$shard_id,
    replications, directory
  )
  run_batch_shard(
    args, configuration, environment, "reference", args$shard_id,
    replications, directory
  )
  message("Shard ", args$shard_id, " completed and validated.")
}

subset_checkpoint <- function(checkpoint, keep) {
  list(
    diagnostics = checkpoint$diagnostics[keep, , drop = FALSE],
    integrals = checkpoint$integrals[keep, , drop = FALSE]
  )
}

combine_checkpoints <- function(checkpoints) {
  diagnostics <- do.call(rbind, lapply(checkpoints, `[[`, "diagnostics"))
  integrals <- do.call(rbind, lapply(checkpoints, `[[`, "integrals"))
  ordering <- order(diagnostics$n, diagnostics$replication)
  diagnostics <- diagnostics[ordering, , drop = FALSE]
  integrals <- integrals[ordering, , drop = FALSE]
  rownames(diagnostics) <- NULL
  rownames(integrals) <- NULL
  list(diagnostics = diagnostics, integrals = integrals)
}

canonical_checkpoint_paths <- function(args, configuration) {
  prefix <- paste0(configuration$run_name, "_")
  list(
    evaluation = file.path(
      args$output_dir, paste0(prefix, "checkpoint.rds")
    ),
    reference = file.path(
      args$output_dir,
      paste0(prefix, "independent_reference_checkpoint.rds")
    ),
    evaluation_metadata = file.path(
      args$output_dir, paste0(prefix, "metadata.rds")
    ),
    reference_metadata = file.path(
      args$output_dir,
      paste0(prefix, "independent_reference_metadata.rds")
    ),
    manifest = file.path(
      args$output_dir, paste0(prefix, "hpc_merge_manifest.rds")
    )
  )
}

validate_base_metadata <- function(metadata, args, configuration, batch,
                                   path) {
  equal_value <- function(x, y) {
    isTRUE(all.equal(x, y, tolerance = 1e-12))
  }
  minimum_replications <- args$rep_start - 1L
  common <-
    identical(metadata$master_seed, args$seed) &&
    equal_value(metadata$epsilon, configuration$epsilon) &&
    equal_value(metadata$rho, configuration$rho) &&
    identical(metadata$cost, "fixed squared Euclidean distance") &&
    identical(metadata$normalized_plan, FALSE) &&
    equal_value(metadata$marginals, configuration$marginals) &&
    equal_value(metadata$population$mass, configuration$population$mass) &&
    equal_value(metadata$population$mean, configuration$population$mean) &&
    equal_value(
      metadata$population$covariance,
      configuration$population$covariance
    ) &&
    !is.null(metadata$target_replications) &&
    metadata$target_replications >= minimum_replications

  if (batch == "evaluation") {
    compatible <-
      common &&
      identical(metadata$n_grid, configuration$n_grid) &&
      identical(
        metadata$test_functions$feature_names,
        configuration$test_functions$feature_names
      ) &&
      equal_value(
        metadata$test_functions$population_integrals,
        configuration$test_functions$population_integrals
      )
  } else {
    compatible <-
      common &&
      identical(metadata$batch, "independent_large_n_reference") &&
      identical(metadata$n, configuration$reference_n) &&
      identical(metadata$seed_stream, 1L) &&
      identical(
        metadata$seed_formula_version,
        "cell_seed_v2_stream_offset_1000000007"
      ) &&
      identical(
        metadata$feature_names,
        configuration$test_functions$feature_names
      ) &&
      equal_value(
        metadata$population_integrals,
        configuration$test_functions$population_integrals
      )
  }
  if (!isTRUE(compatible)) {
    stop(
      "The canonical ", batch, " metadata are incompatible: ", path,
      call. = FALSE
    )
  }
  invisible(TRUE)
}

load_base_checkpoint <- function(args, configuration, batch, path,
                                 metadata_path) {
  if (args$rep_start == 1L) {
    if (file.exists(path)) {
      stop(
        "Canonical checkpoint already exists: ", path,
        ". Use a fresh --output-dir; merge mode never overwrites an existing ",
        "rep-start=1 run.",
        call. = FALSE
      )
    }
    return(empty_checkpoint(configuration$test_functions$feature_names))
  }
  if (!file.exists(path)) {
    stop(
      "--rep-start=", args$rep_start,
      " requires the canonical base checkpoint: ", path, call. = FALSE
    )
  }
  if (!file.exists(metadata_path)) {
    stop(
      "--rep-start=", args$rep_start,
      " requires the companion metadata file: ", metadata_path,
      call. = FALSE
    )
  }
  metadata <- readRDS(metadata_path)
  validate_base_metadata(
    metadata, args, configuration, batch, metadata_path
  )
  checkpoint <- readRDS(path)
  seed_stream <- if (batch == "evaluation") 0L else 1L
  allowed_n <- if (batch == "evaluation") {
    configuration$n_grid
  } else {
    configuration$reference_n
  }
  validate_checkpoint(
    checkpoint, configuration$test_functions$feature_names,
    allowed_n, args$seed, seed_stream, paste0("Canonical ", batch)
  )
  keep <- checkpoint$diagnostics$replication < args$rep_start
  base <- subset_checkpoint(checkpoint, keep)
  base_replications <- seq_len(args$rep_start - 1L)
  expected_design <- make_batch_design(
    configuration, base_replications, batch
  )
  manifest <- make_shard_manifest(
    args, configuration, batch, 0L, base_replications
  )
  manifest$shard_id <- 0L
  manifest$num_shards <- args$num_shards
  object <- list(
    manifest = manifest,
    environment = list(signature = "canonical_base"),
    checkpoint = base
  )
  # Validate the base with the same numerical and design checks, but its
  # manifest is synthetic and used only inside this process.
  validate_complete_batch(
    object, manifest, configuration, expected_design,
    paste0("Canonical ", batch, " base")
  )
  base
}

validate_merged_checkpoint <- function(checkpoint, args, configuration,
                                       batch) {
  seed_stream <- if (batch == "evaluation") 0L else 1L
  allowed_n <- if (batch == "evaluation") {
    configuration$n_grid
  } else {
    configuration$reference_n
  }
  validate_checkpoint(
    checkpoint, configuration$test_functions$feature_names,
    allowed_n, args$seed, seed_stream, paste0("Merged ", batch)
  )
  expected_design <- make_batch_design(
    configuration, seq_len(args$reps), batch
  )
  expected_keys <- job_keys(expected_design)
  observed_keys <- job_keys(checkpoint$diagnostics)
  if (length(observed_keys) != length(expected_keys) ||
      !setequal(observed_keys, expected_keys)) {
    stop("Merged ", batch, " does not contain the exact final design.",
         call. = FALSE)
  }
  if (!all(checkpoint$diagnostics$converged) ||
      any(checkpoint$diagnostics$final_error > 1e-9) ||
      !all(is.finite(checkpoint$integrals)) ||
      any(checkpoint$diagnostics$plan_mass <= 0)) {
    stop("Merged ", batch, " contains an invalid solver result.",
         call. = FALSE)
  }
  invisible(TRUE)
}

run_merge_mode <- function(args, configuration) {
  require_packages(c("MASS", "ggplot2", "reticulate"))
  directory <- configuration_directory(args, configuration)
  if (!dir.exists(directory)) {
    stop("Shard directory does not exist: ", directory, call. = FALSE)
  }
  paths <- canonical_checkpoint_paths(args, configuration)
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  base_evaluation <- load_base_checkpoint(
    args, configuration, "evaluation", paths$evaluation,
    paths$evaluation_metadata
  )
  base_reference <- load_base_checkpoint(
    args, configuration, "reference", paths$reference,
    paths$reference_metadata
  )

  evaluation_parts <- list(base_evaluation)
  reference_parts <- list(base_reference)
  input_paths <- character()
  environment_signatures <- list()
  for (shard_id in seq_len(args$num_shards)) {
    replications <- assigned_replications(args, shard_id)
    for (batch in c("evaluation", "reference")) {
      path <- shard_file(directory, batch, shard_id, args$num_shards)
      if (!file.exists(path)) {
        stop("Missing ", batch, " shard: ", path, call. = FALSE)
      }
      object <- readRDS(path)
      manifest <- make_shard_manifest(
        args, configuration, batch, shard_id, replications
      )
      design <- make_batch_design(configuration, replications, batch)
      label <- paste0(batch, " shard ", shard_id, "/", args$num_shards)
      validate_complete_batch(object, manifest, configuration, design, label)
      if (batch == "evaluation") {
        evaluation_parts[[length(evaluation_parts) + 1L]] <- object$checkpoint
      } else {
        reference_parts[[length(reference_parts) + 1L]] <- object$checkpoint
      }
      input_paths <- c(input_paths, path)
      environment_signatures[[length(environment_signatures) + 1L]] <-
        object$environment$signature
    }
  }
  reference_signature <- environment_signatures[[1L]]
  if (!all(vapply(
    environment_signatures,
    function(x) isTRUE(all.equal(x, reference_signature, tolerance = 0)),
    logical(1)
  ))) {
    stop("Shard R/Python/package version signatures are not identical.",
         call. = FALSE)
  }

  evaluation <- combine_checkpoints(evaluation_parts)
  reference <- combine_checkpoints(reference_parts)
  validate_merged_checkpoint(evaluation, args, configuration, "evaluation")
  validate_merged_checkpoint(reference, args, configuration, "reference")
  if (anyDuplicated(evaluation$diagnostics$seed) ||
      anyDuplicated(reference$diagnostics$seed) ||
      length(intersect(
        evaluation$diagnostics$seed, reference$diagnostics$seed
      ))) {
    stop("Merged evaluation/reference seeds are duplicated or overlap.",
         call. = FALSE)
  }

  atomic_save_rds(evaluation, paths$evaluation)
  atomic_save_rds(reference, paths$reference)
  merge_manifest <- list(
    schema_version = "clt_hpc_merge_v1",
    run_name = configuration$run_name,
    target_replications = args$reps,
    rep_start = args$rep_start,
    num_shards = args$num_shards,
    master_seed = args$seed,
    evaluation_rows = nrow(evaluation$diagnostics),
    reference_rows = nrow(reference$diagnostics),
    feature_names = configuration$test_functions$feature_names,
    shard_md5 = tools::md5sum(input_paths),
    environment_signature = reference_signature,
    source_provenance = source_provenance(),
    output_md5 = tools::md5sum(c(paths$evaluation, paths$reference)),
    analysis_status = "pending",
    merged_at = Sys.time()
  )
  atomic_save_rds(merge_manifest, paths$manifest)

  message("Merged and validated canonical checkpoints:")
  message("  ", paths$evaluation)
  message("  ", paths$reference)
  message("Merge manifest: ", paths$manifest)
  message("Regenerating metadata, summaries, covariance, and plots.")
  main()
  merge_manifest$analysis_status <- "complete"
  merge_manifest$analysis_completed_at <- Sys.time()
  merge_manifest$output_md5 <- tools::md5sum(
    c(
      paths$evaluation, paths$reference,
      paths$evaluation_metadata, paths$reference_metadata
    )
  )
  atomic_save_rds(merge_manifest, paths$manifest)
  message("HPC merge and final analysis completed successfully.")
}

probe_writable_directory <- function(path, label) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) {
    stop("Cannot create ", label, " directory: ", path, call. = FALSE)
  }
  temporary <- tempfile(pattern = ".clt-write-probe-", tmpdir = path)
  saveRDS(list(ok = TRUE), temporary)
  if (!file.exists(temporary)) {
    stop("Cannot write to ", label, " directory: ", path, call. = FALSE)
  }
  unlink(temporary)
  invisible(TRUE)
}

run_preflight <- function(args) {
  require_packages(c("MASS", "ggplot2", "reticulate"))
  record <- environment_record()
  validate_worker_environment(record)
  probe_writable_directory(args$shard_dir, "shard")
  probe_writable_directory(args$output_dir, "output")
  probe_writable_directory(args$plot_dir, "plot")
  marginals <- make_marginals()
  population <- gaussian_uot_population(
    marginals$a, marginals$b, marginals$A, marginals$B, 0.5, 1.0
  )
  if (!isTRUE(all.equal(
    population$mass, 0.588149213661207, tolerance = 1e-12
  ))) {
    stop("Gaussian UOT population golden-mass check failed.", call. = FALSE)
  }
  X0 <- matrix(c(-0.5, 0.1, 0.2, -0.3), nrow = 2, byrow = TRUE)
  X1 <- matrix(c(0.4, -0.2, -0.1, 0.5), nrow = 2, byrow = TRUE)
  solution <- compute_uot_plan(X0, X1, 0.5, 1.0)
  if (!all(dim(solution$plan) == c(2L, 2L)) ||
      !all(is.finite(solution$plan)) || sum(solution$plan) <= 0 ||
      !isTRUE(solution$converged) || !is.finite(solution$final_error) ||
      solution$final_error > 1e-9) {
    stop("POT preflight solve failed.", call. = FALSE)
  }
  cat("HPC preflight passed.\n")
  cat("Population mass:", format(population$mass, digits = 15), "\n")
  cat("Python:", record$Python_executable, "\n")
  print(record$signature)
  invisible(record)
}

main_hpc <- function() {
  args <- parse_hpc_arguments(commandArgs(trailingOnly = TRUE))
  if (args$help) {
    hpc_usage()
    return(invisible(NULL))
  }
  validate_hpc_arguments(args)
  assert_repository_root()
  if (args$mode == "preflight") {
    run_preflight(args)
    return(invisible(NULL))
  }
  configuration <- make_hpc_configuration(args)
  if (args$mode == "shard") {
    run_shard_mode(args, configuration)
  } else {
    run_merge_mode(args, configuration)
  }
  invisible(NULL)
}

if (sys.nframe() == 0L) main_hpc()
