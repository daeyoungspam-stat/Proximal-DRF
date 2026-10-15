`%||%` <- function(x, y) if (is.null(x)) y else x

assert_true <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}

project_root <- function() {
  env_root <- Sys.getenv("PDR_NULL_ROOT", unset = "")
  if (nzchar(env_root)) return(normalizePath(env_root, mustWork = TRUE))
  candidates <- unique(c(getwd(), file.path(getwd(), "Parametric null study", "simulation"),
                         dirname(getwd()), dirname(dirname(getwd()))))
  hit <- candidates[file.exists(file.path(candidates, "config", "simulation_config.R"))]
  if (!length(hit)) stop("Cannot locate simulation root; set PDR_NULL_ROOT.", call. = FALSE)
  normalizePath(hit[[1L]], mustWork = TRUE)
}

source_simulation_files <- function(root = project_root(), compile_cpp = TRUE) {
  source(file.path(root, "config", "simulation_config.R"), local = .GlobalEnv)
  files <- c("dgp.R", "bridges.R", "bandwidth.R", "operators.R", "fast_apply.R",
             "wild_bootstrap.R", "checkpoint.R", "simulation.R", "execution.R")
  for (file in files) source(file.path(root, "R", file), local = .GlobalEnv)
  arc_execution <- file.path(root, "R", "arc_execution.R")
  if (file.exists(arc_execution)) source(arc_execution, local = .GlobalEnv)
  if (compile_cpp) {
    load_fast_cv(root)
    load_fold_fits(root)
    load_fast_fold_apply(root)
  }
  invisible(TRUE)
}

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  answer <- list(); i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token, call. = FALSE)
    token <- substring(token, 3L)
    if (grepl("=", token, fixed = TRUE)) {
      pieces <- strsplit(token, "=", fixed = TRUE)[[1L]]
      answer[[pieces[[1L]]]] <- paste(pieces[-1L], collapse = "=")
    } else if (i < length(args) && !startsWith(args[[i + 1L]], "--")) {
      answer[[token]] <- args[[i + 1L]]; i <- i + 1L
    } else answer[[token]] <- TRUE
    i <- i + 1L
  }
  answer
}

cli_value <- function(args, name, default = NULL,
                      mode = c("character", "integer", "numeric", "logical")) {
  mode <- match.arg(mode); value <- args[[name]] %||% default
  if (is.null(value)) return(NULL)
  switch(mode, character = as.character(value), integer = as.integer(value),
         numeric = as.numeric(value),
         logical = if (is.logical(value)) value else
           tolower(value) %in% c("1", "true", "yes", "y"))
}

seed_for <- function(n, rep, stream = 0L, base_seed = 20260825L) {
  modulus <- 2147483646
  value <- (as.double(base_seed) + 104729 * as.double(n) +
              1009 * as.double(rep) + 7919 * as.double(stream)) %% modulus
  as.integer(value + 1)
}

with_preserved_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed); force(expr)
}

balanced_folds <- function(n, k, seed) {
  assert_true(k >= 2L && k <= n, "Invalid fold count")
  with_preserved_seed(seed, {
    fold <- rep(seq_len(k), length.out = n); fold[sample.int(n)] <- fold; fold
  })
}

safe_solve <- function(a, b = NULL, ridge = 0, condition_limit = 1e12) {
  a <- as.matrix(a); if (ridge > 0) diag(a) <- diag(a) + ridge
  condition <- tryCatch(kappa(a, exact = TRUE), error = function(e) Inf)
  if (!is.finite(condition) || condition > condition_limit)
    stop(sprintf("Ill-conditioned matrix (kappa=%.3e)", condition), call. = FALSE)
  if (is.null(b)) solve(a) else solve(a, b)
}

local_inverse <- function(a, condition_limit = 1e12) {
  a <- as.matrix(a)
  condition <- tryCatch(kappa(a, exact = TRUE), error = function(e) Inf)
  if (!is.finite(condition) || condition > condition_limit)
    stop(sprintf("Ill-conditioned local design (kappa=%.3e)", condition), call. = FALSE)
  list(inverse = solve(a), condition = condition)
}

atomic_save_rds <- function(object, path, compress = "xz") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = compress)
  if (is.null(readRDS(temporary))) stop("Atomic-save validation failed", call. = FALSE)
  if (!file.rename(temporary, path)) stop("Atomic rename failed: ", path, call. = FALSE)
  invisible(path)
}

canonicalize_fingerprint_value <- function(x) {
  # Numerical rules assembled by eigensolvers can differ in their final bits
  # across BLAS/LAPACK implementations.  The quadrature rule is identified by
  # its declarative specification below, so derived arrays are not fingerprint
  # inputs.  Other non-integer numerics are serialized at fourteen significant
  # digits to make the remaining configuration text platform stable.
  if (is.double(x)) {
    answer <- ifelse(
      is.na(x), "NA", ifelse(
        is.nan(x), "NaN", ifelse(
          is.infinite(x), ifelse(x > 0, "Inf", "-Inf"),
          formatC(x, digits = 14L, format = "g", flag = "#")
        )
      )
    )
    names(answer) <- names(x)
    return(answer)
  }
  if (is.data.frame(x)) {
    answer <- lapply(x, canonicalize_fingerprint_value)
    class(answer) <- "data.frame"
    attr(answer, "row.names") <- .set_row_names(nrow(x))
    return(answer)
  }
  if (is.list(x)) return(lapply(x, canonicalize_fingerprint_value))
  x
}

fingerprint_config_view <- function(cfg) {
  derived_quadrature <- c(
    "quadrature_node", "quadrature_weight", "test_weight",
    "weighted_quadrature"
  )
  answer <- cfg[setdiff(names(cfg), derived_quadrature)]
  canonicalize_fingerprint_value(answer)
}

fingerprint_source_files <- function(root = project_root()) {
  # Bind simulation, recovery, postprocessing, reporting, verification, and
  # validation to the same immutable submission.  A source edit at any stage
  # therefore cannot be silently combined with artifacts from another build.
  execution_scripts <- list.files(
    file.path(root, "scripts"), "[.]R$", full.names = TRUE)
  validation_entrypoints <- list.files(
    file.path(root, "tests"), "[.](R|sh)$", full.names = TRUE)
  arc_entrypoints <- list.files(
    file.path(root, "arc"), "[.](R|sh|slurm)$", full.names = TRUE)
  sort(c(list.files(file.path(root, "config"), "[.]R$", full.names = TRUE),
         list.files(file.path(root, "R"), "[.]R$", full.names = TRUE),
         list.files(file.path(root, "src"), "[.]cpp$", full.names = TRUE),
         file.path(root, "run_local_12cores.sh"),
         list.files(file.path(root, "reference"), "[.](rds|csv|json|txt)$", full.names = TRUE),
         execution_scripts, validation_entrypoints, arc_entrypoints))
}

stable_text_hash <- function(text) {
  raw <- as.integer(charToRaw(enc2utf8(paste(text, collapse = "\n"))))
  h <- 0
  for (value in raw) h <- (h * 131 + value) %% 2147483647
  sprintf("%08x", as.integer(h))
}

source_fingerprint_manifest <- function(root = project_root()) {
  files <- fingerprint_source_files(root)
  assert_true(all(file.exists(files)), "Fingerprint manifest contains a missing file")
  files <- normalizePath(files, mustWork = TRUE)
  prefix <- paste0(normalizePath(root, mustWork = TRUE), .Platform$file.sep)
  relative <- ifelse(startsWith(files, prefix), substring(files, nchar(prefix) + 1L), files)
  assert_true(all(startsWith(files, prefix)),
              "Fingerprint source is outside the simulation root")
  data.frame(
    path = relative,
    md5 = unname(tools::md5sum(files)),
    bytes = unname(as.numeric(file.info(files)$size)),
    stringsAsFactors = FALSE
  )
}

source_fingerprint <- function(root = project_root()) {
  manifest <- source_fingerprint_manifest(root)
  stable_text_hash(paste(manifest$path, manifest$md5, manifest$bytes, sep = "|"))
}

config_fingerprint <- function(cfg, root = project_root()) {
  source_manifest <- source_fingerprint_manifest(root)
  manifest <- paste(source_manifest$path, source_manifest$md5,
                    source_manifest$bytes, sep = "|", collapse = "\n")
  fingerprint_cfg <- fingerprint_config_view(cfg)
  text <- paste(paste(capture.output(dput(fingerprint_cfg)), collapse = "\n"), manifest,
                sep = "\n--SOURCE--\n")
  stable_text_hash(text)
}

RUN_MANIFEST_SCHEMA_VERSION <- "pdr-null-immutable-run-manifest-v1"

atomic_save_rds_exclusive <- function(object, path, compress = "xz") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lock <- paste0(path, ".creation_lock")
  if (!dir.create(lock, recursive = FALSE, showWarnings = FALSE))
    stop("Run-manifest creation lock already exists: ", lock, call. = FALSE)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (file.exists(path)) stop("Immutable file already exists: ", path, call. = FALSE)
  temporary <- tempfile(
    pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = compress)
  if (is.null(readRDS(temporary))) stop("Run-manifest validation failed", call. = FALSE)
  if (file.exists(path)) stop("Immutable file appeared during creation: ", path,
                              call. = FALSE)
  if (!file.rename(temporary, path))
    stop("Atomic run-manifest rename failed: ", path, call. = FALSE)
  invisible(path)
}

create_run_manifest <- function(path, run_id, profile, cfg,
                                root = project_root()) {
  assert_true(length(run_id) == 1L && grepl("^[A-Za-z0-9_.-]+$", run_id),
              "Invalid run_id")
  assert_true(identical(as.character(cfg$profile), as.character(profile)),
              "Run-manifest profile does not match the configuration")
  sources <- source_fingerprint_manifest(root)
  answer <- list(
    schema_version = RUN_MANIFEST_SCHEMA_VERSION,
    run_id = as.character(run_id),
    profile = as.character(profile),
    config_id = config_fingerprint(cfg, root),
    source_fingerprint = source_fingerprint(root),
    source_manifest = sources,
    fingerprint_config = fingerprint_config_view(cfg),
    expected_replications = length(cfg$sample_sizes) * cfg$repetitions,
    expected_bootstrap_draws = cfg$bootstrap_draws,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  atomic_save_rds_exclusive(answer, path)
  answer
}

assert_run_manifest <- function(path, run_id, expected_config_id, profile, cfg,
                                root = project_root(), stage = "unspecified") {
  if (!file.exists(path))
    stop("Missing immutable run manifest at stage ", stage, ": ", path,
         call. = FALSE)
  answer <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(answer, "error") || !is.list(answer))
    stop("Unreadable immutable run manifest at stage ", stage, call. = FALSE)
  required <- c(
    "schema_version", "run_id", "profile", "config_id",
    "source_fingerprint", "source_manifest", "fingerprint_config",
    "expected_replications", "expected_bootstrap_draws", "created_utc"
  )
  missing <- setdiff(required, names(answer))
  if (length(missing))
    stop("Run manifest lacks fields: ", paste(missing, collapse = ", "),
         call. = FALSE)
  current_config_id <- config_fingerprint(cfg, root)
  current_source_fingerprint <- source_fingerprint(root)
  current_sources <- source_fingerprint_manifest(root)
  checks <- c(
    schema = identical(answer$schema_version, RUN_MANIFEST_SCHEMA_VERSION),
    run_id = identical(as.character(answer$run_id), as.character(run_id)),
    profile = identical(as.character(answer$profile), as.character(profile)),
    environment_config_id = identical(as.character(answer$config_id),
                                      as.character(expected_config_id)),
    current_config_id = identical(as.character(answer$config_id), current_config_id),
    current_source_fingerprint = identical(
      as.character(answer$source_fingerprint), current_source_fingerprint),
    source_manifest = identical(answer$source_manifest, current_sources),
    config_view = identical(answer$fingerprint_config,
                            fingerprint_config_view(cfg)),
    expected_replications = identical(
      as.integer(answer$expected_replications),
      as.integer(length(cfg$sample_sizes) * cfg$repetitions)),
    bootstrap_draws = identical(
      as.integer(answer$expected_bootstrap_draws), as.integer(cfg$bootstrap_draws))
  )
  if (!all(checks))
    stop("Immutable run-manifest mismatch at stage ", stage, ": ",
         paste(names(checks)[!checks], collapse = ", "), call. = FALSE)
  answer
}

resolve_run_manifest_from_env <- function(root = project_root()) {
  relative <- Sys.getenv("PDR_RUN_MANIFEST_REL", unset = "")
  if (!nzchar(relative)) stop("PDR_RUN_MANIFEST_REL is not set", call. = FALSE)
  components <- strsplit(relative, "[/\\\\]+")[[1L]]
  if (startsWith(relative, "/") || any(components %in% c("", ".", "..")))
    stop("PDR_RUN_MANIFEST_REL must be a normalized relative path", call. = FALSE)
  root <- normalizePath(root, mustWork = TRUE)
  path <- normalizePath(file.path(root, relative), mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if (!startsWith(path, prefix))
    stop("Run manifest resolves outside the simulation root", call. = FALSE)
  path
}

assert_run_manifest_from_env <- function(profile, cfg, root = project_root(),
                                         stage = "unspecified") {
  run_id <- Sys.getenv("PDR_RUN_ID", unset = "")
  expected_config_id <- Sys.getenv("PDR_RUN_CONFIG_ID", unset = "")
  if (!nzchar(run_id) || !nzchar(expected_config_id))
    stop("PDR_RUN_ID and PDR_RUN_CONFIG_ID must be set", call. = FALSE)
  assert_run_manifest(
    resolve_run_manifest_from_env(root), run_id, expected_config_id,
    profile, cfg, root, stage)
}

weighted_integral <- function(value, cfg) sum(cfg$weighted_quadrature * value)

trapz_grid <- function(x, y) {
  ord <- order(x); x <- x[ord]; y <- y[ord]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}
