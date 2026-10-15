# A checkpoint contains completed condition/scenario columns, never a changed
# seed, a changed estimator, or a replacement Monte Carlo replication.
PDR_CHECKPOINT_SCHEMA <- "pdr-condition-checkpoint-v1"

checkpoint_error <- function(message, path = "") {
  stop(structure(list(message = message, call = NULL, path = path),
    class = c("pdr_checkpoint_error", "error", "condition")))
}

checkpoint_pause <- function(context, reason) {
  stop(structure(list(message = reason, call = NULL,
    checkpoint_path = context$path, completed_columns = context$completed,
    warnings = context$warnings),
    class = c("pdr_checkpoint_pause", "error", "condition")))
}

checkpoint_digest <- function(object, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  path <- tempfile(".checkpoint_digest_", tmpdir = directory)
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, compress = FALSE, version = 3L)
  unname(tools::md5sum(path))
}

checkpoint_identity <- function(n, rep, cfg, config_id) list(
  schema = PDR_CHECKPOINT_SCHEMA, config_id = config_id,
  n = as.integer(n), rep = as.integer(rep),
  dgp_seed = seed_for(n, rep, 1L, cfg$base_seed),
  fold_seed = seed_for(n, rep, 2L, cfg$base_seed),
  multiplier_seed = seed_for(n, rep, 500L, cfg$base_seed),
  bootstrap_draws = cfg$bootstrap_draws,
  bootstrap_schemes = cfg$bootstrap_schemes,
  multiplier_laws = cfg$multiplier_laws)

checkpoint_read <- function(path, identity) {
  if (!file.exists(path)) return(NULL)
  envelope <- tryCatch(readRDS(path), error = function(e)
    checkpoint_error(paste("Unreadable checkpoint:", conditionMessage(e)), path))
  if (!is.list(envelope) || !is.list(envelope$payload) ||
      !identical(envelope$schema, PDR_CHECKPOINT_SCHEMA) ||
      !identical(envelope$md5, checkpoint_digest(envelope$payload, dirname(path))))
    checkpoint_error("Checkpoint checksum or envelope is invalid; retained.", path)
  x <- envelope$payload
  if (!identical(x$identity, identity))
    checkpoint_error("Checkpoint identity differs from requested run; retained.", path)
  required <- c("identity", "key", "prep_signature", "completed", "results",
    "bandwidth_rows", "residual_rows", "warnings", "timing", "segments")
  if (!all(required %in% names(x)) || !is.data.frame(x$key) ||
      !all(c("condition", "scenario") %in% names(x$key)) ||
      anyNA(x$key) || anyDuplicated(x$key) ||
      !is.character(x$prep_signature) || length(x$prep_signature) != 1L ||
      is.na(x$prep_signature) || !grepl("^[a-f0-9]{32}$", x$prep_signature) ||
      !is.integer(x$completed) || length(x$completed) != 1L ||
      is.na(x$completed) || x$completed < 0L || x$completed > nrow(x$key) ||
      !is.character(x$warnings) || anyNA(x$warnings) ||
      !is.numeric(x$timing) ||
      !identical(names(x$timing), c("nuisance", "bandwidth", "inference", "total")) ||
      any(!is.finite(x$timing)) || any(x$timing < 0) ||
      !is.integer(x$segments) || length(x$segments) != 1L ||
      is.na(x$segments) || x$segments < 1L)
    checkpoint_error("Checkpoint metadata/structure is invalid; retained.", path)
  for (field in c("results", "bandwidth_rows", "residual_rows")) {
    value <- x[[field]]
    if (!is.list(value) || length(value) != nrow(x$key) ||
        (x$completed > 0L && !all(vapply(value[seq_len(x$completed)],
                                        is.data.frame, logical(1L)))) ||
        (x$completed < nrow(x$key) && !all(vapply(
          value[seq.int(x$completed + 1L, nrow(x$key))], is.null, logical(1L)))))
      checkpoint_error(paste("Invalid completed prefix in", field, "; retained."), path)
  }
  x
}

checkpoint_context <- function(n, rep, cfg, root) {
  e <- new.env(parent = emptyenv())
  e$path <- Sys.getenv("PDR_CHECKPOINT_FILE", unset = "")
  e$enabled <- nzchar(e$path)
  e$config_id <- config_fingerprint(cfg, root)
  expected <- Sys.getenv("PDR_CHECKPOINT_CONFIG_ID", unset = "")
  if (e$enabled && nzchar(expected) && !identical(expected, e$config_id))
    checkpoint_error("Checkpoint configuration changed before replication.", e$path)
  e$identity <- checkpoint_identity(n, rep, cfg, e$config_id)
  e$saved <- if (e$enabled) checkpoint_read(e$path, e$identity) else NULL
  e$completed <- if (is.null(e$saved)) 0L else e$saved$completed
  e$warnings <- if (is.null(e$saved)) character() else e$saved$warnings
  e$prior_timing <- if (is.null(e$saved))
    c(nuisance = 0, bandwidth = 0, inference = 0, total = 0) else e$saved$timing
  e$segments <- if (is.null(e$saved)) 1L else e$saved$segments + 1L
  e
}

checkpoint_validate_columns <- function(saved, cfg, path) {
  template <- result_cell_template(cfg)
  for (j in seq_len(saved$completed)) {
    key <- saved$key[j, , drop = FALSE]
    x <- saved$results[[j]]
    target <- template[template$condition == key$condition &
                         template$scenario == key$scenario, , drop = FALSE]
    numeric_fields <- c("statistic", "bootstrap_critical", "p_value", "bandwidth")
    if (!all(c(names(target), numeric_fields, "reject") %in% names(x)) ||
        anyDuplicated(test_cell_key(x)) ||
        !identical(sort(test_cell_key(x)), sort(test_cell_key(target))) ||
        any(!is.finite(as.matrix(x[numeric_fields]))) ||
        anyNA(x$reject) || any(x$reject != (x$p_value <= cfg$alpha)) ||
        any(x$p_value < 1 / (cfg$bootstrap_draws + 1) | x$p_value > 1) ||
        any(x$bandwidth <= 0))
      checkpoint_error(sprintf("Invalid checkpoint test column %d; retained.", j), path)
    target <- target[match(test_cell_key(x), test_cell_key(target)), , drop = FALSE]
    for (field in names(target)) if (!identical(x[[field]], target[[field]]))
      checkpoint_error(sprintf("Checkpoint column %d field %s differs; retained.",
                               j, field), path)
    band <- saved$bandwidth_rows[[j]]; residual <- saved$residual_rows[[j]]
    if (!all(c("condition", "scenario", "method", "bandwidth") %in% names(band)) ||
        !identical(band$method, cfg$bandwidth_methods) ||
        any(band$condition != key$condition | band$scenario != key$scenario) ||
        any(!is.finite(band$bandwidth) | band$bandwidth <= 0) ||
        !all(c("condition", "scenario") %in% names(residual)) ||
        nrow(residual) != 1L || residual$condition != key$condition ||
        residual$scenario != key$scenario)
      checkpoint_error(sprintf("Invalid checkpoint diagnostics column %d; retained.", j), path)
  }
  invisible(TRUE)
}

checkpoint_prepare <- function(context, key, signature_data, cfg) {
  if (!context$enabled) return(invisible(NULL))
  context$key <- key
  context$prep_signature <- checkpoint_digest(signature_data, dirname(context$path))
  if (!is.null(context$saved)) {
    if (!identical(key, context$saved$key) ||
        !identical(context$prep_signature, context$saved$prep_signature))
      checkpoint_error("Regenerated inputs differ from checkpoint; retained.", context$path)
    checkpoint_validate_columns(context$saved, cfg, context$path)
  }
  invisible(TRUE)
}

checkpoint_commit <- function(context, column, results, bandwidth_rows,
                              residual_rows, current_timing) {
  if (!context$enabled) return(invisible(NULL))
  payload <- list(identity = context$identity, key = context$key,
    prep_signature = context$prep_signature, completed = as.integer(column),
    results = results, bandwidth_rows = bandwidth_rows, residual_rows = residual_rows,
    warnings = unique(context$warnings),
    timing = context$prior_timing + current_timing, segments = context$segments)
  envelope <- list(schema = PDR_CHECKPOINT_SCHEMA, payload = payload,
    md5 = checkpoint_digest(payload, dirname(context$path)))
  atomic_save_rds(envelope, context$path, compress = FALSE)
  context$completed <- as.integer(column)
  invisible(TRUE)
}

checkpoint_check_deadline <- function(context) {
  if (!context$enabled) return(invisible(FALSE))
  raw <- Sys.getenv("PDR_DEADLINE_UNIX", unset = "")
  if (nzchar(raw)) {
    deadline <- suppressWarnings(as.numeric(raw))
    if (length(deadline) != 1L || !is.finite(deadline) || deadline <= 0)
      checkpoint_error("Invalid PDR_DEADLINE_UNIX.", context$path)
    if (as.numeric(Sys.time()) >= deadline)
      checkpoint_pause(context, "Soft time limit reached; completed columns retained.")
  }
  # Validation-only hook: never consulted unless both flags are explicitly set.
  if (identical(Sys.getenv("PDR_CHECKPOINT_TEST_MODE"), "1")) {
    after <- suppressWarnings(as.integer(Sys.getenv("PDR_CHECKPOINT_TEST_AFTER", "")))
    if (length(after) == 1L && !is.na(after) && after >= 0L &&
        context$completed >= after)
      checkpoint_pause(context, "Test-only condition checkpoint pause.")
  }
  invisible(FALSE)
}
