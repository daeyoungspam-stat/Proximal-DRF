RETRY_POLICY_VERSION <- "post-array-fixed-seed-index-v1"
RETRY_MAX_REPLACEMENT_ATTEMPTS <- 5L

retry_seed_index <- function(logical_rep, attempt, repetitions) {
  logical_rep <- as.integer(logical_rep)
  attempt <- as.integer(attempt)
  repetitions <- as.integer(repetitions)
  assert_true(length(logical_rep) == 1L && logical_rep >= 1L &&
                logical_rep <= repetitions, "Invalid logical replication slot")
  assert_true(length(attempt) == 1L && attempt >= 0L,
              "Invalid retry attempt")
  logical_rep + attempt * repetitions
}

retry_seed_triplet <- function(n, logical_rep, attempt, cfg) {
  seed_index <- retry_seed_index(logical_rep, attempt, cfg$repetitions)
  c(
    dgp_seed = seed_for(n, seed_index, 1L, cfg$base_seed),
    fold_seed = seed_for(n, seed_index, 2L, cfg$base_seed),
    multiplier_seed = seed_for(n, seed_index, 500L, cfg$base_seed)
  )
}

assert_retry_seed_uniqueness <- function(cfg,
                                         max_attempt = RETRY_MAX_REPLACEMENT_ATTEMPTS) {
  rows <- expand.grid(
    n = cfg$sample_sizes,
    logical_rep = seq_len(cfg$repetitions),
    attempt = 0:as.integer(max_attempt),
    stream = c(1L, 2L, 500L),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  rows$seed_index <- rows$logical_rep + rows$attempt * cfg$repetitions
  rows$seed <- mapply(
    seed_for, rows$n, rows$seed_index, rows$stream,
    MoreArgs = list(base_seed = cfg$base_seed)
  )
  assert_true(!anyDuplicated(rows[c("n", "seed")]),
              "Retry seed map contains a collision within a sample size")
  invisible(rows)
}

retry_run_dir <- function(outdir, config_id) {
  file.path(outdir, paste0("config_", config_id))
}

retry_primary_failure_path <- function(canonical_path) {
  sub("[.]rds$", ".failed.rds", canonical_path)
}

retry_attempt_dir <- function(outdir, config_id, attempt) {
  file.path(retry_run_dir(outdir, config_id), "arc_retry", "attempts",
            sprintf("attempt_%02d", as.integer(attempt)))
}

retry_attempt_failure_path <- function(outdir, config_id, n, logical_rep, attempt) {
  file.path(retry_attempt_dir(outdir, config_id, attempt), paste0("n_", n),
            sprintf("rep_%04d.failed.rds", as.integer(logical_rep)))
}

retry_candidate_path <- function(outdir, config_id, n, logical_rep, attempt) {
  file.path(retry_attempt_dir(outdir, config_id, attempt), "candidates",
            paste0("n_", n), sprintf("rep_%04d.rds", as.integer(logical_rep)))
}

read_rds_status <- function(path) {
  if (!file.exists(path)) return(list(exists = FALSE, readable = FALSE, value = NULL,
                                      error = ""))
  value <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(value, "error")) return(list(
    exists = TRUE, readable = FALSE, value = NULL, error = conditionMessage(value)))
  list(exists = TRUE, readable = TRUE, value = value, error = "")
}

failure_attempt <- function(failure, default = 0L) {
  metadata <- failure$metadata %||% list()
  as.integer(metadata$retry_attempt %||% metadata$attempt %||% default)
}

failure_classes <- function(failure) {
  metadata <- failure$metadata %||% list()
  unique(as.character(c(
    failure$failure_class %||% character(),
    failure$condition_classes %||% character(),
    failure$error_classes %||% character(),
    metadata$failure_class %||% character(),
    metadata$condition_classes %||% character(),
    metadata$error_classes %||% character()
  )))
}

classify_retry_failure <- function(failure) {
  if (!is.list(failure)) return("unreadable_failure")
  classes <- failure_classes(failure)
  if (any(classes %in% c("pdr_sample_admissibility_error",
                         "sample_admissibility_error"))) {
    return("sample_admissibility")
  }
  # Backward-compatible recognition of the one audited legacy condition.  New
  # code should use the typed classes above; no other message is seed-retried.
  error <- as.character(failure$error %||% "")
  if (identical(error, "Too few observations in local test window"))
    return("sample_admissibility")
  "unexpected_failure"
}

empty_failure_ledger <- function() data.frame(
  n = integer(), rep = integer(), path = character(),
  failure_path = character(), artifact_present = logical(),
  expected_artifact = logical(), valid_part = logical(),
  mapping_source = character(), resolution = character(), unresolved = logical(),
  artifact_config_id = character(), seed_index = integer(),
  retry_attempt = integer(), failure_class = character(), error = character(),
  stringsAsFactors = FALSE
)

failure_artifact_ledger <- function(failure_artifacts, expected, config_id) {
  required <- c("n", "rep", "path", "valid")
  missing <- setdiff(required, names(expected))
  if (length(missing))
    stop("Expected-slot table lacks: ", paste(missing, collapse = ", "),
         call. = FALSE)
  if (!length(failure_artifacts)) return(empty_failure_ledger())
  expected_key <- paste(as.integer(expected$n), as.integer(expected$rep), sep = "|")
  if (anyDuplicated(expected_key))
    stop("Expected replication slots are not unique", call. = FALSE)

  scalar_integer <- function(x) {
    value <- suppressWarnings(as.numeric(x %||% NA_real_))
    if (length(value) != 1L || !is.finite(value) || value != as.integer(value))
      return(NA_integer_)
    as.integer(value)
  }
  scalar_character <- function(x) {
    value <- as.character(x %||% NA_character_)
    if (length(value) != 1L) NA_character_ else value
  }

  rows <- lapply(sort(failure_artifacts), function(failure_path) {
    failure <- tryCatch(readRDS(failure_path), error = function(e) e)
    readable <- !inherits(failure, "error") && is.list(failure)
    read_error <- if (inherits(failure, "error")) conditionMessage(failure) else
      if (!is.list(failure)) "Failure artifact is not an R list" else ""
    metadata <- if (readable) failure$metadata %||% list() else list()
    artifact_config_id <- scalar_character(metadata$config_id)
    metadata_n <- scalar_integer(metadata$n)
    metadata_rep <- scalar_integer(metadata$logical_rep %||% metadata$rep)
    metadata_key <- paste(metadata_n, metadata_rep, sep = "|")
    metadata_mappable <- readable && !is.na(artifact_config_id) &&
      identical(artifact_config_id, as.character(config_id)) &&
      !is.na(metadata_n) && !is.na(metadata_rep)
    expected_index <- if (metadata_mappable)
      match(metadata_key, expected_key) else NA_integer_
    mapped <- !is.na(expected_index)
    valid_part <- mapped && isTRUE(expected$valid[[expected_index]])
    resolution <- if (!mapped) "orphan_failure_artifact" else if (valid_part)
      "resolved_by_valid_part" else if (!readable)
        "unresolved_unreadable_artifact" else "unresolved_failure_artifact"
    classes <- if (readable) failure_classes(failure) else character()
    data.frame(
      n = if (mapped) as.integer(expected$n[[expected_index]]) else metadata_n,
      rep = if (mapped) as.integer(expected$rep[[expected_index]]) else metadata_rep,
      path = if (mapped) expected$path[[expected_index]] else "",
      failure_path = failure_path, artifact_present = TRUE,
      expected_artifact = mapped, valid_part = valid_part,
      mapping_source = if (mapped) "metadata" else "none",
      resolution = resolution, unresolved = resolution != "resolved_by_valid_part",
      artifact_config_id = artifact_config_id,
      seed_index = scalar_integer(metadata$seed_index),
      retry_attempt = scalar_integer(
        metadata$retry_attempt %||% metadata$attempt %||% 0L),
      failure_class = if (length(classes)) classes[[1L]] else "",
      error = if (!readable) read_error else scalar_character(failure$error %||% ""),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

slot_failure_artifacts <- function(outdir, config_id, n, logical_rep,
                                   max_attempt = RETRY_MAX_REPLACEMENT_ATTEMPTS) {
  canonical <- result_path(outdir, config_id, n, logical_rep)
  paths <- c(retry_primary_failure_path(canonical), vapply(
    0:as.integer(max_attempt), function(attempt)
      retry_attempt_failure_path(outdir, config_id, n, logical_rep, attempt),
    character(1L)
  ))
  defaults <- c(0L, 0:as.integer(max_attempt))
  present <- file.exists(paths)
  if (!any(present)) return(data.frame())
  rows <- lapply(which(present), function(j) {
    status <- read_rds_status(paths[[j]])
    failure <- status$value
    attempt <- if (status$readable) failure_attempt(failure, defaults[[j]]) else defaults[[j]]
    data.frame(
      path = paths[[j]], attempt = attempt, readable = status$readable,
      classification = if (status$readable) classify_retry_failure(failure) else
        "unreadable_failure",
      error = if (status$readable) as.character(failure$error %||% "") else status$error,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

read_primary_arc_status <- function(outdir, config_id) {
  status_dir <- file.path(retry_run_dir(outdir, config_id), "arc_status")
  paths <- sort(list.files(
    status_dir, pattern = "^task_[0-9]+[.]csv$", full.names = TRUE))
  if (!length(paths)) return(data.frame())
  rows <- lapply(paths, function(path) {
    value <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(value) || !all(c("n", "rep", "status", "detail") %in% names(value)))
      return(NULL)
    value$task_status_path <- path
    value
  })
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (!length(rows)) return(data.frame())
  answer <- do.call(rbind, rows)
  key <- paste(answer$n, answer$rep, sep = "|")
  answer[!duplicated(key, fromLast = TRUE), , drop = FALSE]
}

primary_status_row <- function(primary_status, n, logical_rep) {
  if (is.null(primary_status) || !nrow(primary_status)) return(NULL)
  hit <- which(primary_status$n == n & primary_status$rep == logical_rep)
  if (!length(hit)) return(NULL)
  primary_status[hit[[length(hit)]], , drop = FALSE]
}

assess_retry_slot <- function(outdir, config_id, n, logical_rep, cfg,
                              max_attempt = RETRY_MAX_REPLACEMENT_ATTEMPTS,
                              primary_status = NULL) {
  canonical <- result_path(outdir, config_id, n, logical_rep)
  if (validate_part(canonical, cfg, n, logical_rep)) return(data.frame(
    n = n, rep = logical_rep, state = "valid", action = "none", attempt = NA_integer_,
    seed_index = NA_integer_, detail = "", stringsAsFactors = FALSE))

  canonical_status <- read_rds_status(canonical)
  failures <- slot_failure_artifacts(
    outdir, config_id, n, logical_rep, max_attempt = max_attempt)
  if (nrow(failures)) {
    failures <- failures[order(failures$attempt, failures$path), , drop = FALSE]
    latest <- failures[nrow(failures), , drop = FALSE]
    if (latest$classification == "sample_admissibility") {
      next_attempt <- latest$attempt + 1L
      if (next_attempt <= max_attempt) return(data.frame(
        n = n, rep = logical_rep, state = "sample_admissibility_failure",
        action = "replacement_seed", attempt = next_attempt,
        seed_index = retry_seed_index(logical_rep, next_attempt, cfg$repetitions),
        detail = latest$error, stringsAsFactors = FALSE))
      return(data.frame(
        n = n, rep = logical_rep, state = "retry_cap_exhausted", action = "hard_stop",
        attempt = next_attempt, seed_index = NA_integer_, detail = latest$error,
        stringsAsFactors = FALSE))
    }
    return(data.frame(
      n = n, rep = logical_rep, state = latest$classification,
      action = "hard_stop", attempt = latest$attempt,
      seed_index = NA_integer_, detail = latest$error, stringsAsFactors = FALSE))
  }

  primary <- primary_status_row(primary_status, n, logical_rep)
  if (!is.null(primary)) {
    primary_state <- as.character(primary$status[[1L]])
    primary_detail <- as.character(primary$detail[[1L]] %||% "")
    if (primary_state == "failed" &&
        identical(primary_detail, "Too few observations in local test window")) {
      return(data.frame(
        n = n, rep = logical_rep, state = "sample_admissibility_status_only",
        action = "replacement_seed", attempt = 1L,
        seed_index = retry_seed_index(logical_rep, 1L, cfg$repetitions),
        detail = primary_detail, stringsAsFactors = FALSE))
    }
    if (primary_state %in% c(
      "runner_error", "source_changed", "invalid_after_save", "failed")) {
      return(data.frame(
        n = n, rep = logical_rep, state = paste0("primary_", primary_state),
        action = "hard_stop", attempt = 0L, seed_index = NA_integer_,
        detail = primary_detail, stringsAsFactors = FALSE))
    }
  }

  if (!canonical_status$exists || !canonical_status$readable) return(data.frame(
    n = n, rep = logical_rep,
    state = if (canonical_status$exists) "unreadable_canonical" else "missing_infrastructure",
    action = "same_seed", attempt = 0L,
    seed_index = retry_seed_index(logical_rep, 0L, cfg$repetitions),
    detail = canonical_status$error, stringsAsFactors = FALSE))

  metadata <- canonical_status$value$metadata %||% list()
  if (!identical(as.character(metadata$config_id %||% ""), config_id)) {
    detail <- "Readable canonical artifact has a different config_id"
  } else {
    detail <- "Readable canonical artifact failed schema/integrity validation"
  }
  data.frame(
    n = n, rep = logical_rep, state = "invalid_canonical", action = "hard_stop",
    attempt = NA_integer_, seed_index = NA_integer_, detail = detail,
    stringsAsFactors = FALSE)
}

assess_all_retry_slots <- function(outdir, config_id, cfg,
                                   max_attempt = RETRY_MAX_REPLACEMENT_ATTEMPTS,
                                   primary_status = NULL) {
  grid <- expand.grid(
    n = cfg$sample_sizes, rep = seq_len(cfg$repetitions),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  rows <- lapply(seq_len(nrow(grid)), function(j) assess_retry_slot(
    outdir, config_id, grid$n[[j]], grid$rep[[j]], cfg,
    max_attempt = max_attempt, primary_status = primary_status))
  do.call(rbind, rows)
}

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."),
                        tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  write.csv(x, temporary, row.names = FALSE)
  check <- tryCatch(read.csv(temporary, stringsAsFactors = FALSE),
                    error = function(e) NULL)
  if (is.null(check) || nrow(check) != nrow(x))
    stop("Atomic CSV validation failed", call. = FALSE)
  if (file.exists(path)) stop("Immutable ledger already exists: ", path, call. = FALSE)
  if (!file.rename(temporary, path)) stop("Atomic CSV rename failed: ", path, call. = FALSE)
  invisible(path)
}

condition_audit_fields <- function(error) {
  fields <- setdiff(names(error), c("message", "call"))
  context <- if (length(fields)) lapply(fields, function(field) error[[field]]) else list()
  names(context) <- fields
  list(condition_classes = class(error), condition_context = context)
}

rewrite_logical_replication <- function(answer, logical_rep, seed_index, attempt,
                                        actual_seeds) {
  answer$metadata$rep <- as.integer(logical_rep)
  answer$metadata$logical_rep <- as.integer(logical_rep)
  answer$metadata$seed_index <- as.integer(seed_index)
  answer$metadata$retry_attempt <- as.integer(attempt)
  answer$metadata$retry_policy_version <- RETRY_POLICY_VERSION
  answer$metadata$dgp_seed <- unname(as.integer(actual_seeds[["dgp_seed"]]))
  answer$metadata$fold_seed <- unname(as.integer(actual_seeds[["fold_seed"]]))
  answer$metadata$multiplier_seed <- unname(as.integer(actual_seeds[["multiplier_seed"]]))
  components <- c("test", "bandwidth", "residual_diagnostics",
                  "nuisance_diagnostics", "dgp_diagnostics", "timing")
  for (component in components) {
    if (is.data.frame(answer[[component]]) && "rep" %in% names(answer[[component]]))
      answer[[component]]$rep <- as.integer(logical_rep)
  }
  answer
}

retry_status_row <- function(n, logical_rep, attempt, seed_index, status,
                             seconds, path, detail = "", failure_class = "",
                             seeds = c(dgp_seed = NA_integer_, fold_seed = NA_integer_,
                                       multiplier_seed = NA_integer_)) {
  data.frame(
    n = as.integer(n), rep = as.integer(logical_rep), attempt = as.integer(attempt),
    seed_index = as.integer(seed_index), dgp_seed = as.integer(seeds[["dgp_seed"]]),
    fold_seed = as.integer(seeds[["fold_seed"]]),
    multiplier_seed = as.integer(seeds[["multiplier_seed"]]),
    status = status, seconds = as.numeric(seconds), path = path,
    failure_class = failure_class, detail = detail,
    stringsAsFactors = FALSE
  )
}

execute_retry_attempt <- function(n, logical_rep, attempt, cfg, root, outdir,
                                  config_id, worker = "retry") {
  canonical <- result_path(outdir, config_id, n, logical_rep)
  seed_index <- retry_seed_index(logical_rep, attempt, cfg$repetitions)
  seeds <- retry_seed_triplet(n, logical_rep, attempt, cfg)
  if (validate_part(canonical, cfg, n, logical_rep)) return(retry_status_row(
    n, logical_rep, attempt, seed_index, "skipped_valid", 0, canonical,
    seeds = seeds))

  dir.create(dirname(canonical), recursive = TRUE, showWarnings = FALSE)
  lock <- acquire_lock(canonical)
  if (is.null(lock)) return(retry_status_row(
    n, logical_rep, attempt, seed_index, "locked", 0, canonical,
    detail = "Canonical artifact lock exists", seeds = seeds))
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (validate_part(canonical, cfg, n, logical_rep)) return(retry_status_row(
    n, logical_rep, attempt, seed_index, "skipped_valid", 0, canonical,
    seeds = seeds))
  if (!identical(config_fingerprint(cfg, root), config_id)) return(retry_status_row(
    n, logical_rep, attempt, seed_index, "source_changed", 0, canonical,
    detail = "Source/configuration changed before retry", seeds = seeds))

  warnings_seen <- character()
  start <- proc.time()[[3L]]
  answer <- tryCatch(withCallingHandlers(
    run_one_replication(n, seed_index, cfg, root),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  elapsed <- proc.time()[[3L]] - start
  failure_path <- retry_attempt_failure_path(
    outdir, config_id, n, logical_rep, attempt)
  if (inherits(answer, "error")) {
    audit <- condition_audit_fields(answer)
    failure <- list(
      metadata = list(
        config_id = config_id, n = as.integer(n), rep = as.integer(logical_rep),
        logical_rep = as.integer(logical_rep), retry_attempt = as.integer(attempt),
        seed_index = as.integer(seed_index), worker = worker,
        dgp_seed = seeds[["dgp_seed"]], fold_seed = seeds[["fold_seed"]],
        multiplier_seed = seeds[["multiplier_seed"]],
        retry_policy_version = RETRY_POLICY_VERSION
      ),
      failure_class = audit$condition_classes[[1L]],
      condition_classes = audit$condition_classes,
      condition_context = audit$condition_context,
      error = conditionMessage(answer), warnings = unique(warnings_seen),
      session = capture.output(sessionInfo())
    )
    atomic_save_rds(failure, failure_path)
    return(retry_status_row(
      n, logical_rep, attempt, seed_index, "failed", elapsed, canonical,
      detail = conditionMessage(answer), failure_class = failure$failure_class,
      seeds = seeds))
  }
  if (!identical(config_fingerprint(cfg, root), config_id)) return(retry_status_row(
    n, logical_rep, attempt, seed_index, "source_changed", elapsed, canonical,
    detail = "Source/configuration changed during retry", seeds = seeds))

  answer <- rewrite_logical_replication(
    answer, logical_rep, seed_index, attempt, seeds)
  answer$warnings <- unique(warnings_seen)
  candidate <- retry_candidate_path(outdir, config_id, n, logical_rep, attempt)
  atomic_save_rds(answer, candidate)
  if (!validate_part(candidate, cfg, n, logical_rep)) {
    failure <- list(
      metadata = list(
        config_id = config_id, n = as.integer(n), rep = as.integer(logical_rep),
        logical_rep = as.integer(logical_rep), retry_attempt = as.integer(attempt),
        seed_index = as.integer(seed_index), worker = worker,
        dgp_seed = seeds[["dgp_seed"]], fold_seed = seeds[["fold_seed"]],
        multiplier_seed = seeds[["multiplier_seed"]],
        retry_policy_version = RETRY_POLICY_VERSION
      ),
      failure_class = "invalid_retry_candidate",
      condition_classes = "invalid_retry_candidate",
      error = "Retry candidate failed schema/integrity validation",
      warnings = unique(warnings_seen), session = capture.output(sessionInfo())
    )
    atomic_save_rds(failure, failure_path)
    return(retry_status_row(
      n, logical_rep, attempt, seed_index, "invalid_after_save", elapsed, candidate,
      detail = failure$error, failure_class = failure$failure_class, seeds = seeds))
  }

  quarantined <- ""
  if (file.exists(canonical)) {
    quarantined <- paste0(canonical, sprintf(".invalid_before_retry_%02d", attempt))
    if (file.exists(quarantined) || !file.rename(canonical, quarantined))
      return(retry_status_row(
        n, logical_rep, attempt, seed_index, "promotion_failed", elapsed, candidate,
        detail = "Could not quarantine invalid canonical artifact", seeds = seeds))
  }
  if (!file.rename(candidate, canonical)) return(retry_status_row(
    n, logical_rep, attempt, seed_index, "promotion_failed", elapsed, candidate,
    detail = "Could not promote validated retry candidate", seeds = seeds))
  if (!validate_part(canonical, cfg, n, logical_rep)) return(retry_status_row(
    n, logical_rep, attempt, seed_index, "invalid_after_promotion", elapsed, canonical,
    detail = "Promoted canonical artifact failed validation", seeds = seeds))
  retry_status_row(
    n, logical_rep, attempt, seed_index, "completed", elapsed, canonical,
    detail = if (nzchar(quarantined)) paste0("quarantined=", quarantined) else "",
    seeds = seeds)
}
