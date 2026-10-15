args<-commandArgs(TRUE); request<-readRDS(args[[1L]]); destination<-args[[2L]]
options(warn=1)
digest<-function(x) { f<-tempfile(); on.exit(unlink(f)); saveRDS(x,f,version=3L); unname(tools::md5sum(f)) }
if(isTRUE(request$prepare)) {
  engines<-c(if(any(c("joint","proxlq") %in% request$components)) "joint",
    if("null" %in% request$components) c("null_original","null_recovery5","null_recovery35","null_recovery13"))
  for(engine in engines) {
    input<-tempfile(fileext=".rds")
    subrequest<-request;subrequest$prepare<-FALSE;subrequest$compile_engine<-engine
    saveRDS(subrequest,input)
    code<-system2(file.path(R.home("bin"),"Rscript"),c("--vanilla",
      shQuote(file.path(request$assets,"runner.R")),shQuote(input),shQuote(destination)))
    unlink(input)
    if(code!=0L) stop("Compilation failed for ",engine)
  }
  saveRDS(engines,destination);quit(status=0L)
}
if(!is.null(request$compile_engine)) {
  root<-file.path(request$assets,"engines",request$compile_engine)
  cache<-file.path(request$cache,request$compile_engine)
  dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  Sys.setenv(PDR_JOINT_ROOT=root,PDR_NULL_ROOT=root,
    PDR_JOINT_CPP_CACHE=cache,PDR_PROXLQ_CPP_CACHE=cache)
  source(file.path(root,"R/utils.R"));source_simulation_files(root,TRUE)
  if(request$compile_engine=="joint" && "proxlq" %in% request$components) {
    source(file.path(root,"proxlq_addon/R/addon.R"));source_proxlq_addon(root,TRUE)
  }
  quit(status=0L)
}
dir.create(request$cache,recursive=TRUE,showWarnings=FALSE)
Sys.setenv(PDR_JOINT_CPP_CACHE=request$cache,PDR_PROXLQ_CPP_CACHE=request$cache)
if(isTRUE(request$validation)) {
  roots<-file.path(request$assets,"engines",c("joint","null_original"))
  scripts<-c(file.path(roots[1],"scripts/validate_math.R"),
    file.path(roots[1],"proxlq_addon/scripts/validate_math.R"),
    file.path(roots[2],"tests/bootstrap_equivalence.R"),
    file.path(roots[2],"tests/fast_apply_validation.R"),
    file.path(request$assets,"engines/null_recovery13/tests/test_cv_admissibility.R"))
  Sys.setenv(PDR_JOINT_ROOT=roots[1],PDR_NULL_ROOT=roots[2])
  codes<-vapply(scripts,function(s) {
    previous<-getwd(); on.exit(setwd(previous))
    setwd(dirname(s))
    if(grepl("null_recovery13",s,fixed=TRUE))
      Sys.setenv(PDR_NULL_ROOT=file.path(request$assets,"engines/null_recovery13"))
    system2(file.path(R.home("bin"),"Rscript"),
      c("--vanilla",shQuote(basename(s)),"--quick"))
  },integer(1))
  if(any(codes!=0L)) stop("Engine validation failed.")
  saveRDS(list(checks=basename(scripts),exit_code=codes),destination)
  quit(status=0L)
}
task<-request$task
engine<-if(task$component=="null") "null_original" else "joint"
if(task$component=="null" && task$profile=="production") {
  if(task$n==5000L && task$rep %in% c(373L,374L,382L,383L,613L)) engine<-"null_recovery5"
  if(task$n==2500L && task$rep %in% c(17L,39L,80L,117L,124L,158L,160L,161L,
    205L,211L,233L,267L,283L,289L,323L,344L,388L,416L,438L,440L,469L,476L,
    482L,509L,588L,617L,619L,621L,623L,640L,677L,700L,709L,750L,788L)) engine<-"null_recovery35"
  if(task$n==2500L && task$rep %in% c(818L,819L,821L,835L,845L,889L,907L,910L,
    918L,925L,935L,983L,986L)) engine<-"null_recovery13"
}
root<-file.path(request$assets,"engines",engine)
cache<-file.path(request$cache,engine)
dir.create(cache,recursive=TRUE,showWarnings=FALSE)
Sys.setenv(PDR_JOINT_CPP_CACHE=cache,PDR_PROXLQ_CPP_CACHE=cache)
Sys.setenv(PDR_JOINT_ROOT=root,PDR_NULL_ROOT=root)
source(file.path(root,"R/utils.R"))
source_simulation_files(root,compile_cpp=TRUE)
cfg<-simulation_config(task$profile)
if(task$component=="proxlq") {
  source(file.path(root,"proxlq_addon/R/addon.R"))
  source_proxlq_addon(root,TRUE)
  cfg<-proxlq_addon_config(task$profile,root)
}
valid<-function(result,path) {
  if(task$component=="null") return(validate_result_object(result,cfg,task$n,task$rep,config_fingerprint(cfg,root)))
  f<-tempfile(fileext=".rds"); on.exit(unlink(f)); saveRDS(result,f)
  if(task$component=="joint") validate_part(f,cfg,task$n,task$rep) else
    validate_proxlq_part(f,cfg,task$n,task$rep,root)
}
if(file.exists(destination)) {
  existing<-readRDS(destination)
  stopifnot(identical(existing$manifest_id,request$manifest_id),
    identical(existing$key,task$key),identical(existing$payload_md5,digest(existing$result)),
    valid(existing$result,destination))
  cat("Verified existing artifact\n"); quit(status=0L)
}
if(task$component=="null") Sys.setenv(PDR_CHECKPOINT_FILE=request$checkpoint)
answer<-if(task$component=="proxlq") run_one_proxlq_replication(task$n,task$rep,cfg,root) else
  run_one_replication(task$n,task$rep,cfg,root)
stopifnot(valid(answer,destination))
envelope<-list(key=task$key,manifest_id=request$manifest_id,engine=engine,
  configuration=cfg,result=answer,payload_md5=digest(answer))
temporary<-paste0(destination,".tmp"); saveRDS(envelope,temporary)
stopifnot(identical(readRDS(temporary)$payload_md5,digest(answer)),file.rename(temporary,destination))
cat("Completed ",task$key,"\n",sep="")
