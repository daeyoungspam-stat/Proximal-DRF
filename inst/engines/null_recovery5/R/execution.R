result_path <- function(outdir, config_id, n, rep) {
  file.path(outdir, paste0("config_", config_id), paste0("n_", n),
            sprintf("rep_%04d.rds", rep))
}

result_cell_template <- function(cfg) {
  conditions <- simulation_conditions(cfg)
  do.call(rbind, lapply(seq_len(nrow(conditions)), function(j)
    do.call(rbind, lapply(cfg$bootstrap_schemes, function(scheme) {
      specs <- bootstrap_scheme_specs(conditions$degree[[j]], scheme, cfg)
      grid <- expand.grid(scenario = cfg$scenarios$scenario,
                          method = specs$method,
                          multiplier_law = cfg$multiplier_laws,
                          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      grid$condition <- conditions$condition[[j]]
      grid$null <- conditions$null[[j]]
      grid$degree <- conditions$degree[[j]]
      grid$delta <- conditions$delta[[j]]
      grid$bootstrap_scheme <- scheme
      grid$identified <- cfg$scenarios$identified[
        match(grid$scenario, cfg$scenarios$scenario)]
      grid
    }))))
}

test_cell_key <- function(x, paired = FALSE) {
  fields <- c("condition", "scenario", "method", "multiplier_law")
  if (!paired) fields <- c(fields, "bootstrap_scheme")
  do.call(paste, c(x[fields], sep = "|"))
}

.validation_cache <- new.env(parent = emptyenv())
validation_template <- function(cfg, config_id) {
  if (!exists(config_id, .validation_cache, inherits = FALSE)) {
    conditions <- simulation_conditions(cfg)
    test <- result_cell_template(cfg)
    band <- expand.grid(condition = conditions$condition,
                        scenario = cfg$scenarios$scenario,
                        method = cfg$bandwidth_methods,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    residual <- unique(band[c("condition", "scenario")])
    nuisance <- expand.grid(fold = seq_len(cfg$folds),
      nuisance = c("h_correct", "h_omit_w", "q_quadratic", "q_main_effects"),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    assign(config_id, list(test = test, test_key = test_cell_key(test),
      band_key = do.call(paste, c(band, sep = "|")),
      residual_key = do.call(paste, c(residual, sep = "|")),
      nuisance_key = do.call(paste, c(nuisance, sep = "|"))), .validation_cache)
  }
  get(config_id, .validation_cache, inherits = FALSE)
}

validate_result_object <- function(answer, cfg, n, rep, config_id) {
  tryCatch({
    if (!is.list(answer) || is.null(answer$metadata)) return(FALSE)
    md <- answer$metadata
    if (!identical(as.character(md$config_id), as.character(config_id)) ||
        !identical(as.integer(md$n), as.integer(n)) ||
        !identical(as.integer(md$rep), as.integer(rep)) ||
        !identical(as.integer(md$bootstrap_draws), as.integer(cfg$bootstrap_draws)) ||
        !identical(as.character(md$bootstrap_schemes), cfg$bootstrap_schemes) ||
        !identical(as.character(md$multiplier_laws), cfg$multiplier_laws)) return(FALSE)
    template <- validation_template(cfg, config_id)
    exact_key <- function(observed, expected)
      length(observed) == length(expected) && !anyNA(observed) &&
      !anyDuplicated(observed) && identical(sort(observed), sort(expected))
    x <- answer$test
    if (!is.data.frame(x) || !exact_key(test_cell_key(x), template$test_key))
      return(FALSE)
    ix <- match(test_cell_key(x), template$test_key)
    if (!all(x$n == n & x$rep == rep) ||
        !all(x$null == template$test$null[ix]) ||
        !all(x$degree == template$test$degree[ix]) ||
        !all(x$delta == template$test$delta[ix]) ||
        !all(x$identified == template$test$identified[ix])) return(FALSE)
    required_finite <- c("statistic", "bootstrap_critical", "p_value", "bandwidth")
    optional_finite <- c("bootstrap_critical_aligned", "null_fit_identity_error",
      "null_smoothing_bias_max", "null_smoothing_bias_l2", "max_condition",
      "min_local_n", "max_reproduction_error", "residual_sd",
      "residual_abs_q99", "residual_abs_max", "residual_excess_kurtosis")
    if (!all(required_finite %in% names(x))) return(FALSE)
    cols <- c(required_finite, intersect(optional_finite, names(x)))
    if (!all(is.finite(as.matrix(x[cols]))) || anyNA(x$reject) ||
        !all(x$p_value >= 1 / (cfg$bootstrap_draws + 1) & x$p_value <= 1) ||
        !all(x$reject == (x$p_value <= cfg$alpha)) || !all(x$bandwidth > 0))
      return(FALSE)
    if ("null_fit_identity_error" %in% names(x) &&
        max(x$null_fit_identity_error) >= 1e-9) return(FALSE)
    if ("max_reproduction_error" %in% names(x) &&
        max(x$max_reproduction_error) > cfg$local_reproduction_tolerance)
      return(FALSE)
    component_key <- function(data, fields)
      do.call(paste, c(data[fields], sep = "|"))
    exact_key(component_key(answer$bandwidth, c("condition", "scenario", "method")),
              template$band_key) &&
      exact_key(component_key(answer$residual_diagnostics, c("condition", "scenario")),
                template$residual_key) &&
      exact_key(component_key(answer$nuisance_diagnostics, c("fold", "nuisance")),
                template$nuisance_key) &&
      is.data.frame(answer$timing) && nrow(answer$timing) == 1L &&
      is.finite(answer$timing$total_seconds) && answer$timing$total_seconds >= 0
  }, error = function(e) FALSE)
}

validate_part <- function(path, cfg, n, rep, expected_config_id = NULL) {
  if (!file.exists(path)) return(FALSE)
  config_id <- expected_config_id %||% config_fingerprint(cfg)
  answer <- tryCatch(readRDS(path), error = function(e) NULL)
  validate_result_object(answer, cfg, n, rep, config_id)
}

local_hostname <- function() unname(Sys.info()[["nodename"]])
local_pid_alive <- function(pid, hostname = local_hostname()) {
  if (!identical(as.character(hostname), local_hostname()) ||
      length(pid) != 1L || !is.finite(pid) || pid < 1) return(NA)
  isTRUE(tryCatch(tools::pskill(as.integer(pid), signal = 0L),
                  error = function(e) FALSE))
}

lock_owner <- function(lock) {
  path <- file.path(lock, "owner.rds")
  if (!file.exists(path)) return(NULL)
  suppressWarnings(tryCatch(readRDS(path), error = function(e) NULL))
}

acquire_lock <- function(path, stale_hours = NULL) {
  lock <- paste0(path, ".lock")
  if (dir.exists(lock)) {
    owner <- lock_owner(lock)
    age <- as.numeric(difftime(Sys.time(), file.info(lock)$mtime, units = "secs"))
    # Long computations retain locks. Only verified dead local PIDs or empty
    # abandoned creation locks can be reclaimed; archives preserve their record.
    reclaim <- if (!is.null(owner)) identical(local_pid_alive(
      owner$pid, owner$hostname), FALSE) else is.finite(age) && age > 300
    if (reclaim) {
      archive <- paste0(lock, ".stale_", format(Sys.time(), "%Y%m%dT%H%M%S"),
                        "_", Sys.getpid())
      if (!file.rename(lock, archive)) return(NULL)
    }
  }
  if (!dir.create(lock, recursive = FALSE, showWarnings = FALSE)) return(NULL)
  atomic_save_rds(list(pid = Sys.getpid(), hostname = local_hostname(),
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
    file.path(lock, "owner.rds"), compress = FALSE)
  lock
}

execute_replication <- function(n, rep, cfg, root, outdir, force = FALSE,
                                worker = "single", expected_config_id = NULL) {
  config_id <- config_fingerprint(cfg, root)
  if (!is.null(expected_config_id) && !identical(config_id, expected_config_id))
    stop("Source or configuration changed after the worker started.", call. = FALSE)
  path <- result_path(outdir, config_id, n, rep)
  status_row <- function(status, seconds = 0, detail = "") data.frame(
    n = n, rep = rep, worker = worker, pid = Sys.getpid(), status = status,
    seconds = seconds, path = path, detail = detail, stringsAsFactors = FALSE)
  if (validate_part(path, cfg, n, rep, config_id))
    return(status_row("skipped_valid"))
  if (file.exists(path)) return(status_row("invalid_existing", detail =
    "Existing artifact fails validation; retained without overwrite."))
  failed_path <- sub("[.]rds$", ".failed.rds", path)
  if (file.exists(failed_path) && !force) return(status_row("skipped_failed", detail =
    "Recorded failure retained; automatic replacement seeds are disabled."))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lock <- acquire_lock(path)
  if (is.null(lock)) return(status_row("locked", detail = "Artifact lock exists"))
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (validate_part(path, cfg, n, rep, config_id)) return(status_row("skipped_valid"))
  checkpoint_path <- paste0(path, ".checkpoint.rds")
  if (identical(Sys.getenv("PDR_ENABLE_CHECKPOINTS"), "1")) {
    env_names <- c("PDR_CHECKPOINT_FILE", "PDR_CHECKPOINT_CONFIG_ID")
    old_env <- Sys.getenv(env_names, unset = NA_character_)
    on.exit({
      for (j in seq_along(env_names)) {
        if (is.na(old_env[[j]])) Sys.unsetenv(env_names[[j]]) else
          do.call(Sys.setenv, setNames(list(old_env[[j]]), env_names[[j]]))
      }
    }, add = TRUE)
    Sys.setenv(PDR_CHECKPOINT_FILE = checkpoint_path,
               PDR_CHECKPOINT_CONFIG_ID = config_id)
  }
  warnings_seen <- character(); start <- proc.time()[[3L]]
  answer <- tryCatch(withCallingHandlers(run_one_replication(n, rep, cfg, root),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w));
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  elapsed <- proc.time()[[3L]] - start
  if (inherits(answer, "pdr_checkpoint_pause"))
    return(status_row("paused", elapsed, sprintf("%s Completed columns: %d; %s",
      conditionMessage(answer), answer$completed_columns, answer$checkpoint_path)))
  if (!inherits(answer, "error") &&
      !identical(config_fingerprint(cfg, root), config_id))
    return(status_row("source_changed", elapsed,
      "Source or configuration changed during the replication; artifact not committed."))
  if (!inherits(answer, "error") &&
      !validate_result_object(answer, cfg, n, rep, config_id)) {
    atomic_save_rds(answer, sub("[.]rds$", ".invalid.rds", path))
    answer <- simpleError("Replication output failed structural/numerical validation")
  }
  if (inherits(answer, "error")) {
    context_names <- setdiff(names(answer), c("message", "call"))
    failure <- list(metadata = list(config_id = config_id, n = n, rep = rep,
      worker = worker, pid = Sys.getpid(), hostname = local_hostname(),
      logical_rep = rep, seed_index = rep, attempt = 0L, retry_attempt = 0L,
      dgp_seed = seed_for(n, rep, 1L, cfg$base_seed),
      fold_seed = seed_for(n, rep, 2L, cfg$base_seed),
      multiplier_seed = seed_for(n, rep, 500L, cfg$base_seed),
      failed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
      failure_class = class(answer)[[1L]], error_classes = class(answer),
      error_context = as.list(answer[context_names]), error = conditionMessage(answer),
      seconds = elapsed, warnings = unique(warnings_seen),
      session = capture.output(sessionInfo()))
    if (file.exists(failed_path)) failed_path <- sub("[.]rds$",
      paste0(".attempt_", format(Sys.time(), "%Y%m%dT%H%M%S"), ".rds"), failed_path)
    atomic_save_rds(failure, failed_path)
    return(status_row("failed", elapsed, conditionMessage(answer)))
  }
  answer$warnings <- unique(c(answer$warnings, warnings_seen))
  answer$metadata$worker <- worker
  answer$metadata$pid <- Sys.getpid()
  answer$metadata$hostname <- local_hostname()
  atomic_save_rds(answer, path)
  # A fully validated, atomically committed result supersedes its partial
  # checkpoint. Never remove a checkpoint after a pause or a failure.
  if (file.exists(checkpoint_path) && validate_part(path, cfg, n, rep, config_id))
    unlink(checkpoint_path)
  status_row("completed", elapsed)
}

run_task_grid <- function(cfg) {
  # n varies fastest, yielding timing information for every n early in a run.
  expand.grid(n = cfg$sample_sizes, rep = seq_len(cfg$repetitions),
              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
}

audit_run <- function(cfg, outdir, config_id, read_timing = TRUE) {
  tasks <- run_task_grid(cfg)
  tasks$status <- "missing"; tasks$seconds <- NA_real_
  tasks$completed_mtime <- NA_real_
  tasks$failure_class <- tasks$failure_message <- ""
  for (j in seq_len(nrow(tasks))) {
    path <- result_path(outdir, config_id, tasks$n[j], tasks$rep[j])
    if (file.exists(path)) {
      part <- tryCatch(readRDS(path), error = function(e) NULL)
      good <- validate_result_object(part, cfg, tasks$n[j], tasks$rep[j], config_id)
      tasks$status[j] <- if (good) "valid" else "invalid"
      if (good && read_timing) tasks$seconds[j] <- part$timing$total_seconds
      if (good) tasks$completed_mtime[j] <- as.numeric(file.info(path)$mtime)
    } else if (file.exists(sub("[.]rds$", ".failed.rds", path))) {
      failure <- tryCatch(readRDS(sub("[.]rds$", ".failed.rds", path)),
                          error = function(e) NULL)
      tasks$status[j] <- "failed"
      tasks$failure_class[j] <- failure$failure_class %||% "unreadable_failure"
      tasks$failure_message[j] <- failure$error %||% "Unreadable failure artifact"
    } else if (file.exists(paste0(path, ".checkpoint.rds"))) {
      checkpoint_path <- paste0(path, ".checkpoint.rds")
      state <- tryCatch({
        saved <- checkpoint_read(checkpoint_path, checkpoint_identity(
          tasks$n[j], tasks$rep[j], cfg, config_id))
        checkpoint_validate_columns(saved, cfg, checkpoint_path)
        saved
      }, error = function(e) e)
      if (inherits(state, "error")) {
        tasks$status[j] <- "invalid_checkpoint"
        tasks$failure_class[j] <- class(state)[[1L]]
        tasks$failure_message[j] <- conditionMessage(state)
      } else {
        tasks$status[j] <- "checkpointed"
        tasks$seconds[j] <- state$timing[["total"]]
      }
    }
    lock <- paste0(path, ".lock")
    if (dir.exists(lock)) {
      owner <- lock_owner(lock)
      alive <- if (is.null(owner)) NA else local_pid_alive(owner$pid, owner$hostname)
      if (tasks$status[j] %in% c("missing", "checkpointed")) tasks$status[j] <-
        if (isTRUE(alive)) "running" else if (identical(alive, FALSE)) "stale_lock"
        else "unknown_lock"
    }
  }
  tasks
}
