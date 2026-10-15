`%||%` <- function(x, y) if (is.null(x)) y else x

assert_true <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}

project_root <- function() {
  env_root <- Sys.getenv("PDR_JOINT_ROOT", unset = "")
  if (nzchar(env_root)) return(normalizePath(env_root, mustWork = TRUE))
  candidates <- unique(c(
    getwd(),
    file.path(getwd(), "Derivative R study", "joint_simulation"),
    dirname(getwd()), dirname(dirname(getwd()))
  ))
  hit <- candidates[file.exists(file.path(candidates, "config", "simulation_config.R"))]
  if (!length(hit)) stop("Cannot locate joint_simulation root; set PDR_JOINT_ROOT.", call. = FALSE)
  normalizePath(hit[[1L]], mustWork = TRUE)
}

source_simulation_files <- function(root = project_root(), compile_cpp = TRUE) {
  source(file.path(root, "config", "simulation_config.R"), local = .GlobalEnv)
  files <- c("dgp.R", "bridges.R", "bandwidth.R", "inference.R",
             "simulation.R", "execution.R")
  for (file in files) source(file.path(root, "R", file), local = .GlobalEnv)
  if (isTRUE(compile_cpp)) load_fast_cv(root)
  invisible(root)
}

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  answer <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token, call. = FALSE)
    token <- substring(token, 3L)
    if (grepl("=", token, fixed = TRUE)) {
      pieces <- strsplit(token, "=", fixed = TRUE)[[1L]]
      answer[[pieces[[1L]]]] <- paste(pieces[-1L], collapse = "=")
    } else if (i < length(args) && !startsWith(args[[i + 1L]], "--")) {
      answer[[token]] <- args[[i + 1L]]
      i <- i + 1L
    } else {
      answer[[token]] <- TRUE
    }
    i <- i + 1L
  }
  answer
}

cli_value <- function(args, name, default = NULL,
                      mode = c("character", "integer", "numeric", "logical")) {
  mode <- match.arg(mode)
  value <- args[[name]] %||% default
  if (is.null(value)) return(NULL)
  switch(mode,
    character = as.character(value),
    integer = as.integer(value),
    numeric = as.numeric(value),
    logical = if (is.logical(value)) value else tolower(value) %in% c("1", "true", "yes", "y")
  )
}

seed_for <- function(n, rep, stream = 0L, base_seed = 20260823L) {
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
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  force(expr)
}

balanced_folds <- function(n, k, seed) {
  assert_true(k >= 2L && k <= n, "Invalid fold count")
  with_preserved_seed(seed, {
    fold <- rep(seq_len(k), length.out = n)
    fold[sample.int(n)] <- fold
    fold
  })
}

safe_solve <- function(a, b = NULL, ridge = 0, condition_limit = 1e12) {
  a <- as.matrix(a)
  if (ridge > 0) diag(a) <- diag(a) + ridge
  condition <- tryCatch(kappa(a, exact = TRUE), error = function(e) Inf)
  if (!is.finite(condition) || condition > condition_limit)
    stop(sprintf("Ill-conditioned matrix (kappa=%.3e)", condition), call. = FALSE)
  if (is.null(b)) solve(a) else solve(a, b)
}

local_inverse <- function(a, condition_limit = 1e12) {
  a <- as.matrix(a)
  condition <- tryCatch(kappa(a, exact = TRUE), error = function(e) Inf)
  if (!is.finite(condition) || condition > condition_limit)
    stop(sprintf("Ill-conditioned unregularized local design (kappa=%.3e)", condition), call. = FALSE)
  list(inverse = solve(a), condition = condition)
}

trapz_grid <- function(x, y) {
  order <- order(x)
  x <- x[order]
  y <- y[order]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

atomic_save_rds <- function(object, path, compress = "xz") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = compress)
  check <- readRDS(temporary)
  if (is.null(check)) stop("Atomic-save validation failed for ", path, call. = FALSE)
  if (!file.rename(temporary, path)) stop("Atomic rename failed for ", path, call. = FALSE)
  invisible(path)
}

config_fingerprint <- function(cfg, root = tryCatch(project_root(), error = function(e) NULL)) {
  manifest <- ""
  if (!is.null(root)) {
    files <- sort(c(
      list.files(file.path(root, "config"), pattern = "[.]R$", full.names = TRUE),
      list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE),
      list.files(file.path(root, "src"), pattern = "[.]cpp$", full.names = TRUE),
      list.files(file.path(root, "scripts"), pattern = "^run_(task|local_worker)[.]R$", full.names = TRUE)
    ))
    manifest <- paste(substring(files, nchar(root) + 2L),
                      unname(tools::md5sum(files)), collapse = "\n")
  }
  text <- paste(
    paste(capture.output(dput(cfg, control = c("keepNA", "keepInteger"))), collapse = "\n"),
    manifest, sep = "\n--SOURCE--\n"
  )
  raw <- as.integer(charToRaw(enc2utf8(text)))
  h <- 0
  for (value in raw) h <- (h * 131 + value) %% 2147483647
  sprintf("%08x", as.integer(h))
}

psd_correlation <- function(r, tolerance = 1e-10) {
  r <- (r + t(r)) / 2
  eigenvalue <- eigen(r, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalue) >= -tolerance) {
    diag(r) <- 1
    return(list(correlation = r, projected = FALSE, projection_norm = 0,
                min_eigen = min(eigenvalue)))
  }
  eig <- eigen(r, symmetric = TRUE)
  values <- pmax(eig$values, tolerance)
  out <- eig$vectors %*% (values * t(eig$vectors))
  diagonal <- sqrt(pmax(diag(out), tolerance))
  out <- out / outer(diagonal, diagonal)
  diag(out) <- 1
  out <- (out + t(out)) / 2
  list(correlation = out, projected = TRUE,
       projection_norm = sqrt(sum((out - r)^2)), min_eigen = min(eigenvalue))
}

grid_membership <- function(x, grid, tolerance = 1e-12) {
  vapply(x, function(value) any(abs(value - grid) <= tolerance), logical(1))
}
