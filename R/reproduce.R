.assets <- function() {
  p <- system.file(package="pdrReproduce")
  if (!nzchar(p)) stop("Install pdrReproduce before execution.")
  p
}

.digest <- function(x) {
  f <- tempfile(); on.exit(unlink(f))
  saveRDS(x,f,version=3L)
  unname(tools::md5sum(f))
}

.engines <- function() {
  m <- readRDS(file.path(.assets(),"engine_manifest.rds"))
  paths <- file.path(.assets(),"engines",m$engine,m$path)
  if (!all(file.exists(paths)) || !identical(unname(tools::md5sum(paths)),m$md5))
    stop("Engine checksum verification failed.")
  m
}

simulation_plan <- function(study=c("all","inference","parametric_null"),
                            profile=c("production","test"), n=NULL, reps=NULL) {
  study <- match.arg(study); profile <- match.arg(profile)
  components <- switch(study,all=c("joint","proxlq","null"),
    inference=c("joint","proxlq"),parametric_null="null")
  rows <- lapply(components,function(component) {
    sizes <- if (!is.null(n)) n else if(profile=="production")
      c(2500L,5000L,7500L,10000L) else 800L
    iterations <- if(!is.null(reps)) reps else seq_len(if(profile=="test") 2L else
      if(component=="null") 1000L else 5000L)
    if (anyNA(c(sizes,iterations)) || any(c(sizes,iterations)<1) ||
        any(c(sizes,iterations)!=as.integer(c(sizes,iterations))) ||
        anyDuplicated(sizes) || anyDuplicated(iterations)) stop("Invalid n or reps.")
    d <- expand.grid(n=as.integer(sizes),rep=as.integer(iterations),KEEP.OUT.ATTRS=FALSE)
    d$component <- component; d$profile <- profile; d
  })
  answer <- do.call(rbind,rows); rownames(answer)<-NULL
  answer <- answer[order(answer$component,answer$n,answer$rep),]
  answer$key <- sprintf("%s_n%d_rep%04d",answer$component,answer$n,answer$rep)
  rownames(answer)<-NULL
  answer
}

.child <- function(request, output, logfile) {
  input <- tempfile(fileext=".rds"); on.exit(unlink(input))
  saveRDS(request,input)
  system2(file.path(R.home("bin"),"Rscript"),
    c("--vanilla",shQuote(file.path(.assets(),"runner.R")),shQuote(input),shQuote(output)),
    stdout=logfile,stderr=logfile,
    env=c("OMP_NUM_THREADS=1","OPENBLAS_NUM_THREADS=1","MKL_NUM_THREADS=1",
          "VECLIB_MAXIMUM_THREADS=1","BLIS_NUM_THREADS=1"))
}

run_simulations <- function(outdir, study=c("all","inference","parametric_null"),
                            profile=c("production","test"), workers=1L,
                            n=NULL,reps=NULL) {
  study<-match.arg(study); profile<-match.arg(profile)
  if(length(workers)!=1L || is.na(workers) || workers<1 || workers!=as.integer(workers))
    stop("workers must be a positive integer.")
  if (.Platform$OS.type=="windows" && workers>1L)
    stop("Parallel execution requires Linux or macOS; select workers=1 on Windows.")
  plan <- simulation_plan(study,profile,n,reps)
  manifest <- list(schema=1L,plan=plan,engines=.engines(),
    runner_md5=unname(tools::md5sum(file.path(.assets(),"runner.R"))))
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  outdir <- normalizePath(outdir,mustWork=TRUE)
  lock <- file.path(outdir,"controller.lock")
  if(!dir.create(lock,showWarnings=FALSE)) stop("Controller lock exists: ",lock)
  completed_controller <- FALSE
  on.exit(if(completed_controller) unlink(lock,recursive=TRUE),add=TRUE)
  saveRDS(list(pid=Sys.getpid(),started=Sys.time()),file.path(lock,"controller.rds"))
  mf <- file.path(outdir,"manifest.rds")
  if(file.exists(mf)) {
    if(!identical(readRDS(mf),manifest)) stop("Existing run has a different plan or source.")
  } else saveRDS(manifest,mf)
  for(d in c("parts","logs","checkpoints","cache"))
    dir.create(file.path(outdir,d),showWarnings=FALSE)
  prepared <- .child(list(assets=.assets(),prepare=TRUE,components=unique(plan$component),
    cache=file.path(outdir,"cache")),file.path(outdir,"prepared.rds"),
    file.path(outdir,"logs","prepare.log"))
  if(prepared!=0L) stop("Numerical compilation failed; inspect logs/prepare.log.")
  run <- function(i) {
    task<-as.list(plan[i,]); path<-file.path(outdir,"parts",paste0(task$key,".rds"))
    request<-list(task=task,assets=.assets(),manifest_id=.digest(manifest),
      checkpoint=file.path(outdir,"checkpoints",paste0(task$key,".rds")),
      cache=file.path(outdir,"cache"))
    status <- .child(request,path,file.path(outdir,"logs",paste0(task$key,".log")))
    data.frame(key=task$key,exit_code=status)
  }
  result <- if(workers==1L) lapply(seq_len(nrow(plan)),run) else
    parallel::mclapply(seq_len(nrow(plan)),run,mc.cores=as.integer(workers),
      mc.preschedule=FALSE,mc.set.seed=FALSE)
  saveRDS(result,file.path(outdir,"execution_status.rds"))
  completed_controller <- TRUE
  answer<-simulation_status(outdir)
  if(any(answer$status!="valid")) warning("Run contains incomplete tasks; inspect logs and status.")
  answer
}

simulation_status <- function(outdir) {
  manifest<-readRDS(file.path(outdir,"manifest.rds")); id<-.digest(manifest)
  d<-manifest$plan
  d$status<-vapply(d$key,function(key) {
    p<-file.path(outdir,"parts",paste0(key,".rds"))
    if(!file.exists(p)) return("missing")
    x<-tryCatch(readRDS(p),error=function(e)NULL)
    if(is.null(x)||!identical(x$manifest_id,id)||!identical(x$key,key)||
      !identical(x$payload_md5,.digest(x$result))) "invalid" else "valid"
  },character(1))
  d
}

validate_engines <- function(outdir=tempfile("pdr-validation-")) {
  .engines()
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  outdir<-normalizePath(outdir,mustWork=TRUE)
  request<-list(assets=.assets(),validation=TRUE,cache=file.path(outdir,"cache"))
  code<-.child(request,file.path(outdir,"checks.rds"),file.path(outdir,"validation.log"))
  if(code!=0L) stop("Validation failed; inspect ",file.path(outdir,"validation.log"))
  readRDS(file.path(outdir,"checks.rds"))
}
