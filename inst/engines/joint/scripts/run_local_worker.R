root <- normalizePath(Sys.getenv("PDR_JOINT_ROOT", unset = ""), mustWork = TRUE)
source(file.path(root, "R", "utils.R"), local = .GlobalEnv)
source_simulation_files(root, compile_cpp = TRUE)
cli <- parse_cli()
profile <- cli_value(cli, "profile", "production", "character")
cfg <- simulation_config(profile)
worker_id <- cli_value(cli, "worker-id", NULL, "integer")
workers <- cli_value(cli, "workers", NULL, "integer")
n_values <- as.integer(strsplit(cli_value(cli, "n", paste(cfg$sample_sizes, collapse = ","),
                                         "character"), ",", fixed = TRUE)[[1L]])
repetitions <- cli_value(cli, "reps", cfg$repetitions, "integer")
outdir <- normalizePath(cli_value(cli, "outdir", file.path(root, "results"), "character"),
                        mustWork = FALSE)
force <- cli_value(cli, "force", FALSE, "logical")
assert_true(worker_id >= 1L && worker_id <= workers, "Invalid worker id")
cfg$sample_sizes <- n_values; cfg$repetitions <- repetitions
validate_requested_design(cfg)
tasks <- expand.grid(n = n_values, rep = seq_len(repetitions),
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
tasks <- tasks[order(tasks$n, tasks$rep), , drop = FALSE]
owned <- which((seq_len(nrow(tasks)) - 1L) %% workers == worker_id - 1L)
status_dir <- file.path(outdir, paste0("config_", config_fingerprint(cfg, root)), "worker_status")
dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(status_dir, sprintf("worker_%02d.csv", worker_id))
if (file.exists(status_path)) unlink(status_path)
for (position in seq_along(owned)) {
  task <- tasks[owned[[position]], ]
  status <- execute_replication(task$n, task$rep, cfg, root, outdir, force,
                                worker = sprintf("local-%02d", worker_id))
  utils::write.table(status, status_path, sep = ",", row.names = FALSE,
                     col.names = !file.exists(status_path), append = file.exists(status_path),
                     quote = TRUE)
  if (position %% 25L == 0L || position == length(owned)) {
    message(sprintf("worker %02d: %d/%d tasks; last=%s n=%d rep=%d",
                    worker_id, position, length(owned), status$status, task$n, task$rep))
  }
}
