load_fast_fold_apply <- function(root = project_root()) {
  if (exists("blas_dense_apply_cpp", mode = "function", inherits = TRUE))
    return(invisible(TRUE))
  assert_true(requireNamespace("Rcpp", quietly = TRUE), "Package 'Rcpp' is required")
  r_config <- function(name) {
    value <- system2(file.path(R.home("bin"), "R"), c("CMD", "config", name),
                      stdout = TRUE, stderr = TRUE)
    assert_true(is.null(attr(value, "status")) || attr(value, "status") == 0L,
                paste("Cannot read R linker configuration", name, paste(value, collapse = " ")))
    paste(value, collapse = " ")
  }
  # Linux module builds can require the configured Fortran runtime at link
  # time. The macOS R framework supplies its shared BLAS runtime itself;
  # adding its optional external gfortran SDK flags is unnecessary here.
  link_flags <- paste(r_config("BLAS_LIBS"),
                      if (identical(Sys.info()[["sysname"]], "Darwin")) "" else r_config("FLIBS"))
  old <- Sys.getenv("PKG_LIBS", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("PKG_LIBS") else Sys.setenv(PKG_LIBS = old), add = TRUE)
  Sys.setenv(PKG_LIBS = paste(if (is.na(old)) "" else old, link_flags))
  cache <- Sys.getenv("PDR_JOINT_CPP_CACHE", unset = "")
  if (!nzchar(cache)) cache <- file.path(root, ".rcpp_cache")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(file.path(root, "src", "fast_fold_apply.cpp"),
                  rebuild = FALSE, showOutput = FALSE, verbose = FALSE,
                  cacheDir = cache)
  invisible(TRUE)
}

fast_fold_operator_apply <- function(operator, response,
                                      backend = c("matrix", "blas")) {
  backend <- match.arg(backend)
  response <- as.matrix(response)
  if (backend == "matrix") return(as.matrix(operator$S %*% response))
  smoother <- operator$dense_S %||% as.matrix(operator$S)
  blas_dense_apply_cpp(smoother, response)
}

fast_fold_backend_diagnostics <- function() {
  info <- sessionInfo()
  variables <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
  list(implementation = "R-linked BLAS dgemm",
       blas = info$BLAS %||% unname(extSoftVersion()["BLAS"]),
       lapack = info$LAPACK %||% NA_character_,
       R = R.version.string, platform = R.version$platform,
       thread_environment = Sys.getenv(variables, unset = NA_character_))
}
