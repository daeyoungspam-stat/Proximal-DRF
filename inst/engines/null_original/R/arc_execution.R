# ARC orchestration only. Statistical estimators, seeds and bootstrap draws are
# supplied by the unchanged simulation engine and configuration.
ARC_POINTER_SCHEMA <- "pdr-fold-residuals-arc-one-node-pointer-v2"

arc_topology <- function(array_tasks = 36L, nodes = 1L, workers = 40L,
                         concurrent_jobs = 9L) {
  x <- as.integer(c(array_tasks, nodes, workers, concurrent_jobs))
  assert_true(all(is.finite(x)) && all(x >= 1L), "Invalid ARC topology")
  assert_true(nodes == 1L && workers <= 40L &&
                nodes * concurrent_jobs <= 9L &&
                nodes * workers * concurrent_jobs <= 360L &&
                concurrent_jobs <= 9L,
              "ARC topology exceeds the declared node/core/job allocation")
  list(array_tasks = x[1], nodes = x[2], workers = x[3],
       concurrent_jobs = x[4], simultaneous_workers = x[2] * x[3] * x[4])
}

arc_task_plan <- function(cfg, topology = arc_topology()) {
  # Scheduling constants only: these never enter the DGP, estimators or seeds.
  # Same-n shards keep all 40 cores occupied for approximately the same duration.
  # Caps target 20--30 hours at the supplied ARC timings, not a runtime guarantee.
  cap <- c(`2500` = 400L, `5000` = 160L, `7500` = 80L, `10000` = 80L)
  hours <- c(`2500` = 2.52, `5000` = 6.5, `7500` = 10.5, `10000` = 14.6)
  shards <- list()
  for (n in cfg$sample_sizes) {
    count <- unname(cap[as.character(n)])
    if (is.na(count)) count <- 40L  # small validation profiles only
    estimate <- unname(hours[as.character(n)])
    if (is.na(estimate)) estimate <- 0
    for (first in seq.int(1L, cfg$repetitions, by = count)) {
      reps <- seq.int(first, min(first + count - 1L, cfg$repetitions))
      shards[[length(shards) + 1L]] <- data.frame(n = n, rep = reps,
        estimated_shard_hours = ceiling(length(reps) / topology$workers) * estimate)
    }
  }
  assert_true(length(shards) <= topology$array_tasks, "Too few arrays for the task plan")
  # Offer longer shards first; Slurm, not this ordering, controls dispatch times.
  order_shards <- order(-vapply(shards, function(z) z$estimated_shard_hours[1L], numeric(1)),
                        -vapply(shards, function(z) z$n[1L], integer(1)),
                        vapply(shards, function(z) z$rep[1L], integer(1)))
  shards <- shards[order_shards]
  tasks <- do.call(rbind, lapply(seq_along(shards), function(j) {
    z <- shards[[j]]; z$array_index <- as.integer(j); z$node_index <- 0L; z
  }))
  rownames(tasks) <- NULL
  tasks$task_index <- seq_len(nrow(tasks))
  assert_true(!anyDuplicated(tasks[c("n", "rep")]), "Duplicate simulation task")
  expected <- run_task_grid(cfg)
  key <- function(z) paste(z$n, z$rep, sep = ":")
  assert_true(identical(sort(key(tasks)), sort(key(expected))), "Missing simulation tasks")
  tasks
}

arc_job_context <- function() {
  job_id <- Sys.getenv("SLURM_JOB_ID", unset = "")
  array_job <- Sys.getenv("SLURM_ARRAY_JOB_ID", unset = "")
  array_task <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
  job_ref <- if (nzchar(array_job) && nzchar(array_task))
    paste0(array_job, "_", array_task) else job_id
  list(job_id = job_id, job_ref = job_ref, array_job_id = array_job,
       array_task_id = array_task, node_index = Sys.getenv("SLURM_NODEID", ""),
       cluster = Sys.getenv("SLURM_CLUSTER_NAME", ""), hostname = local_hostname(),
       pid = Sys.getpid())
}

arc_require_compute <- function(profile) {
  if (nzchar(Sys.getenv("SLURM_JOB_ID", ""))) return(invisible(TRUE))
  allow_test <- identical(profile, "test") &&
    identical(Sys.getenv("PDR_ARC_ALLOW_LOCAL_TEST", ""), "1")
  assert_true(allow_test, "Run this R entry point inside a Slurm compute allocation")
  invisible(TRUE)
}

arc_terminal_states <- c("BOOT_FAIL", "CANCELLED", "COMPLETED", "DEADLINE",
  "FAILED", "NODE_FAIL", "OUT_OF_MEMORY", "PREEMPTED", "REVOKED", "TIMEOUT")
arc_live_states <- c("PENDING", "RUNNING", "CONFIGURING", "COMPLETING",
  "SUSPENDED", "RESIZING", "REQUEUED", "REQUEUE_FED", "REQUEUE_HOLD",
  "RESV_DEL_HOLD", "SIGNALING", "SPECIAL_EXIT", "STAGE_OUT", "STOPPED")

arc_normalize_slurm_state <- function(state) {
  toupper(sub("[+[:space:]].*$", "", trimws(state)))
}

arc_command_lines <- function(command, args) {
  executable <- Sys.which(command)
  if (!nzchar(executable)) return(NULL)
  ans <- suppressWarnings(tryCatch(system2(executable, args, stdout = TRUE,
    stderr = FALSE, timeout = 15L), error = function(e) NULL))
  if (is.null(ans) || (!is.null(attr(ans, "status")) && attr(ans, "status") != 0L))
    return(NULL)
  ans[nzchar(trimws(ans))]
}

.arc_job_status_cache <- new.env(parent = emptyenv())
arc_slurm_job_state <- function(owner, command = arc_command_lines, cache = TRUE) {
  job_ref <- owner$job_ref %||% owner$job_id %||% ""
  job_id <- owner$job_id %||% job_ref
  if (!grepl("^[0-9]+(_[0-9]+)?$", job_ref)) return("UNKNOWN")
  cluster <- owner$cluster %||% ""
  current_cluster <- Sys.getenv("SLURM_CLUSTER_NAME", "")
  if (nzchar(cluster) && nzchar(current_cluster) && cluster != current_cluster)
    return("UNKNOWN")
  key <- paste(cluster, job_ref, sep = "|")
  now <- as.numeric(Sys.time())
  if (cache && exists(key, .arc_job_status_cache, inherits = FALSE)) {
    entry <- get(key, .arc_job_status_cache)
    if (now - entry$time < 15) return(entry$state)
  }
  parse_state <- function(lines) {
    if (is.null(lines) || !length(lines)) return(character())
    pieces <- strsplit(lines, "|", fixed = TRUE)
    id <- trimws(vapply(pieces, function(x) x[1], character(1)))
    state <- vapply(pieces, function(x) if (length(x) >= 2L) x[2] else "", character(1))
    arc_normalize_slurm_state(state[id %in% c(job_ref, job_id)])
  }
  queued <- parse_state(command("squeue", c("--noheader", "--jobs", shQuote(job_ref),
    "--format", shQuote("%i|%T"))))
  state <- if (length(queued) && any(queued %in% arc_live_states))
    queued[which(queued %in% arc_live_states)[1]] else "UNKNOWN"
  if (state == "UNKNOWN") {
    accounted <- parse_state(command("sacct", c("--noheader", "--allocations",
      "--jobs", shQuote(job_ref), "--format=JobIDRaw,State", "--parsable2")))
    if (length(accounted) && any(accounted %in% arc_live_states))
      state <- accounted[which(accounted %in% arc_live_states)[1]]
    else if (length(accounted) && all(accounted %in% arc_terminal_states))
      state <- accounted[1]
  }
  if (cache) assign(key, list(time = now, state = state), .arc_job_status_cache)
  state
}

arc_lock_liveness <- function(owner, slurm_query = arc_slurm_job_state,
                               local_probe = local_pid_alive) {
  if (is.null(owner)) return("unknown")
  same_host <- identical(as.character(owner$hostname), local_hostname())
  if (same_host) {
    alive <- local_probe(owner$pid, owner$hostname)
    if (isTRUE(alive)) return("local_alive")
    if (identical(alive, FALSE)) return("dead")
  }
  state <- slurm_query(owner)
  if (state %in% arc_terminal_states) return("dead")
  if (state %in% arc_live_states) return("job_alive")
  "unknown"
}

arc_acquire_lock <- function(path, stale_hours = NULL) {
  lock <- paste0(path, ".lock")
  # Every acquisition, including creation of a new lock, is serialized through
  # this guard. A contender must not archive a different owner's newly created
  # lock using a liveness decision made before another contender reclaimed it.
  # An abandoned guard is not reclaimed using elapsed time or an inferred PID;
  # retain it for explicit inspection rather than risking concurrent ownership.
  guard <- paste0(lock, ".reclaim_guard")
  if (!dir.create(guard, recursive = FALSE, showWarnings = FALSE)) return(NULL)
  on.exit(unlink(guard, recursive = TRUE), add = TRUE)
  guard_owner <- arc_job_context()
  guard_owner$started_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  atomic_save_rds(guard_owner, file.path(guard, "owner.rds"), compress = FALSE)
  if (dir.exists(lock)) {
    owner <- lock_owner(lock)
    # Terminal states are checked freshly before this destructive rename; a
    # cached terminal state may predate a Slurm requeue with the same job ID.
    liveness <- arc_lock_liveness(owner, slurm_query = function(x)
      arc_slurm_job_state(x, cache = FALSE))
    if (identical(liveness, "dead")) {
      archive <- paste0(lock, ".stale_", format(Sys.time(), "%Y%m%dT%H%M%S"),
        "_", Sys.getenv("SLURM_JOB_ID", "local"), "_", Sys.getpid())
      if (!file.rename(lock, archive)) return(NULL)
    } else return(NULL)
  }
  if (!dir.create(lock, recursive = FALSE, showWarnings = FALSE)) return(NULL)
  owner <- arc_job_context()
  owner$started_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  atomic_save_rds(owner, file.path(lock, "owner.rds"), compress = FALSE)
  lock
}

# execute_replication resolves this name at call time. No estimator code changes.
acquire_lock <- arc_acquire_lock

arc_initialize <- function(cfg, root, outdir, topology = arc_topology()) {
  config_id <- config_fingerprint(cfg, root)
  run_id <- paste0("arc_", config_id)
  run_dir <- file.path(outdir, paste0("config_", config_id))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(run_dir, "run_manifest.rds")
  if (!file.exists(manifest_path))
    invisible(create_run_manifest(manifest_path, run_id, cfg$profile, cfg, root))
  invisible(assert_run_manifest(manifest_path, run_id, config_id, cfg$profile,
                                cfg, root, "arc-initialize"))
  plan <- arc_task_plan(cfg, topology)
  plan_path <- file.path(run_dir, "arc_task_plan.rds")
  if (!file.exists(plan_path)) atomic_save_rds_exclusive(plan, plan_path)
  assert_true(identical(readRDS(plan_path), plan), "Existing ARC task plan differs")
  pointer <- list(schema = ARC_POINTER_SCHEMA, config_id = config_id,
    run_id = run_id, profile = cfg$profile, topology = topology,
    manifest_relative = file.path(paste0("config_", config_id), "run_manifest.rds"),
    expected_replications = nrow(plan), bootstrap_draws = cfg$bootstrap_draws)
  pointer_path <- file.path(outdir, "ARC_RUN_POINTER.rds")
  if (!file.exists(pointer_path)) atomic_save_rds_exclusive(pointer, pointer_path)
  assert_true(identical(readRDS(pointer_path), pointer),
              "This outdir belongs to a different immutable ARC run; use a new outdir")
  pointer
}

arc_assert_run <- function(cfg, root, outdir, stage) {
  pointer_path <- file.path(outdir, "ARC_RUN_POINTER.rds")
  assert_true(file.exists(pointer_path), "ARC initialization pointer is missing")
  pointer <- readRDS(pointer_path)
  assert_true(identical(pointer$schema, ARC_POINTER_SCHEMA), "Invalid ARC pointer schema")
  assert_true(identical(pointer$config_id, config_fingerprint(cfg, root)),
              "ARC source/configuration changed after initialization")
  assert_true(identical(pointer$profile, cfg$profile), "ARC profile mismatch")
  invisible(assert_run_manifest(file.path(outdir, pointer$manifest_relative),
    pointer$run_id, pointer$config_id, cfg$profile, cfg, root, stage))
  Sys.setenv(PDR_RUN_ID = pointer$run_id, PDR_RUN_CONFIG_ID = pointer$config_id)
  pointer
}

arc_task_status <- function(n, rep, cfg, outdir, config_id) {
  path <- result_path(outdir, config_id, n, rep)
  if (validate_part(path, cfg, n, rep, config_id)) return("valid")
  if (file.exists(path)) return("invalid")
  if (file.exists(sub("[.]rds$", ".failed.rds", path))) return("failed")
  lock <- paste0(path, ".lock")
  if (dir.exists(paste0(lock, ".reclaim_guard"))) return("unknown_lock")
  if (dir.exists(lock)) return(switch(arc_lock_liveness(lock_owner(lock)),
    dead = "stale_lock", local_alive = "running", job_alive = "running_job", "unknown_lock"))
  partial <- paste0(path, ".checkpoint.rds")
  if (!file.exists(partial)) return("missing")
  tryCatch({
    saved <- checkpoint_read(partial, checkpoint_identity(n, rep, cfg, config_id))
    checkpoint_validate_columns(saved, cfg, partial)
    "checkpointed"
  }, error = function(e) "invalid_checkpoint")
}

arc_audit_run <- function(cfg, outdir, config_id) {
  audit <- audit_run(cfg, outdir, config_id)
  guarded <- vapply(seq_len(nrow(audit)), function(j) dir.exists(paste0(
    result_path(outdir, config_id, audit$n[j], audit$rep[j]),
    ".lock.reclaim_guard")), logical(1L))
  locked <- which(audit$status %in% c("running", "stale_lock", "unknown_lock") | guarded)
  for (j in locked) audit$status[j] <- arc_task_status(
    audit$n[j], audit$rep[j], cfg, outdir, config_id)
  audit
}

arc_node_tasks <- function(cfg, topology, array_index, node_index) {
  assert_true(array_index >= 1L && array_index <= topology$array_tasks,
              "Array index is outside the initialized plan")
  assert_true(node_index >= 0L && node_index < topology$nodes,
              "Node index is outside the initialized plan")
  plan <- arc_task_plan(cfg, topology)
  plan[plan$array_index == array_index & plan$node_index == node_index, , drop = FALSE]
}

arc_queue_next <- function(queue_path, tasks) {
  lock <- NULL; start <- Sys.time()
  while (is.null(lock)) {
    lock <- acquire_lock(queue_path)
    if (is.null(lock)) {
      if (as.numeric(difftime(Sys.time(), start, units = "secs")) > 120)
        stop("Local task queue lock could not be verified/recovered in120seconds")
      Sys.sleep(0.025)
    }
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  index <- readRDS(queue_path)
  if (index > nrow(tasks)) return(NULL)
  atomic_save_rds(index + 1L, queue_path, compress = FALSE)
  tasks[index, , drop = FALSE]
}
