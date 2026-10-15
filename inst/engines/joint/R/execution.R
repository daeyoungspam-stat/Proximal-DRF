failure_path <- function(path) sub("[.]rds$", ".failed.rds", path)
lock_path <- function(path) paste0(path, ".lock")

acquire_lock <- function(path, stale_hours = 12) {
  lock <- lock_path(path)
  if (dir.exists(lock)) {
    age <- as.numeric(difftime(Sys.time(), file.info(lock)$mtime, units = "hours"))
    if (is.finite(age) && age > stale_hours) unlink(lock, recursive = TRUE)
  }
  if (!dir.create(lock, recursive = FALSE, showWarnings = FALSE)) return(NULL)
  lock
}

execute_replication <- function(n, rep, cfg, root, outdir,
                                force = FALSE, worker = "single") {
  config_id <- config_fingerprint(cfg, root)
  path <- result_path(outdir, config_id, n, rep)
  if (!force && validate_part(path, cfg, n, rep)) {
    return(data.frame(n = n, rep = rep, worker = worker, status = "skipped_valid",
                      seconds = 0, path = path, detail = "", stringsAsFactors = FALSE))
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lock <- acquire_lock(path)
  if (is.null(lock)) {
    return(data.frame(n = n, rep = rep, worker = worker, status = "locked",
                      seconds = 0, path = path, detail = "artifact lock exists",
                      stringsAsFactors = FALSE))
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  warnings_seen <- character(); started <- proc.time()[[3L]]
  answer <- tryCatch(withCallingHandlers(
    run_one_replication(n, rep, cfg, root),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) e)
  elapsed <- proc.time()[[3L]] - started
  if (inherits(answer, "error")) {
    failure <- list(
      metadata = list(config_id = config_id, n = n, rep = rep, worker = worker,
                      failed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)),
      error = conditionMessage(answer),
      call = paste(deparse(conditionCall(answer)), collapse = " "),
      warnings = unique(warnings_seen), session = capture.output(sessionInfo())
    )
    atomic_save_rds(failure, failure_path(path))
    return(data.frame(n = n, rep = rep, worker = worker, status = "failed",
                      seconds = elapsed, path = path, detail = conditionMessage(answer),
                      stringsAsFactors = FALSE))
  }
  answer$warnings <- unique(warnings_seen)
  atomic_save_rds(answer, path)
  if (file.exists(failure_path(path))) unlink(failure_path(path))
  if (!validate_part(path, cfg, n, rep)) {
    return(data.frame(n = n, rep = rep, worker = worker, status = "invalid_after_save",
                      seconds = elapsed, path = path,
                      detail = "saved artifact failed mathematical/schema validation",
                      stringsAsFactors = FALSE))
  }
  data.frame(n = n, rep = rep, worker = worker, status = "completed",
             seconds = elapsed, path = path, detail = "", stringsAsFactors = FALSE)
}
