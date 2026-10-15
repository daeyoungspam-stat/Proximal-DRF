summarize_simulations <- function(outdir) {
  audit<-simulation_status(outdir)
  if(any(audit$status!="valid")) stop("Complete valid task plan required for summaries.")
  accum<-new.env(parent=emptyenv())
  add<-function(group,keys,values) {
    key<-do.call(paste,c(keys,sep="|")); ix<-order(key)
    keys<-keys[ix,,drop=FALSE]; values<-as.matrix(values[ix,,drop=FALSE])
    if(anyDuplicated(key)) stop("Duplicate summary cell within replication.")
    if(any(!is.finite(values))) stop("Nonfinite summary input.")
    if(exists(group,accum,inherits=FALSE)) {
      z<-get(group,accum); stopifnot(identical(z$keys,keys))
      z$sum<-z$sum+values; z$count<-z$count+1L
    } else z<-list(keys=keys,sum=values,count=1L)
    assign(group,z,accum)
  }
  for(i in seq_len(nrow(audit))) {
    x<-readRDS(file.path(outdir,"parts",paste0(audit$key[i],".rds")))$result
    if(audit$component[i]=="null") {
      d<-x$test
      d<-d[d$bootstrap_scheme=="fold_parametric",]
      rownames(d)<-NULL
      add(paste0("test_",audit$n[i]),
        d[c("n","scenario","null","degree","delta","method","multiplier_law","bootstrap_scheme")],
        data.frame(rejection_probability=as.numeric(d$reject),mean_statistic=d$statistic,
          mean_critical=d$bootstrap_critical,mean_bandwidth=d$bandwidth))
    } else {
      d<-x$curve
      if(audit$component[i]=="joint") d<-d[!(d$estimand=="derivative" & d$method=="ProxLL-US"),]
      else {d<-d[d$method=="ProxLQ-US-log10",];d$method<-"ProxLQ-US"}
      rownames(d)<-NULL
      add(paste0("curve_",audit$component[i],"_",audit$n[i]),
        d[c("n","estimand","scenario","method","a","on_point_grid","truth")],
        data.frame(mean_estimate=d$estimate,mean_estimate_squared=d$estimate^2,
          mean_se=d$se,coverage=as.numeric(d$pointwise_coverage),mean_width=d$ci_width,
          mse=d$squared_error,mean_bandwidth=d$bandwidth))
      if(audit$component[i]=="joint") {
        b<-x$band;rownames(b)<-NULL
        add(paste0("band_",audit$n[i]),b[c("n","estimand","scenario","method")],
          data.frame(coverage=as.numeric(b$simultaneous_coverage),mean_width=b$mean_width))
      }
    }
  }
  finish<-function(prefix) {
    groups<-ls(accum,pattern=paste0("^",prefix))
    if(!length(groups)) return(data.frame())
    do.call(rbind,lapply(groups,function(g) {
      z<-get(g,accum); a<-cbind(z$keys,replications=z$count,as.data.frame(z$sum/z$count))
      if(prefix=="curve_") {
        a$bias<-a$mean_estimate-a$truth
        a$empirical_sd<-if(z$count>1L) sqrt(pmax(0,(a$mean_estimate_squared-a$mean_estimate^2)*z$count/(z$count-1))) else NA_real_
        a$rmse<-sqrt(a$mse)
      }
      if(prefix=="test_") a$monte_carlo_se<-sqrt(a$rejection_probability*(1-a$rejection_probability)/z$count)
      a
    }))
  }
  answer<-list(curve=finish("curve_"),band=finish("band_"),test=finish("test_"),audit=audit)
  saveRDS(answer,file.path(outdir,"summary.rds"))
  answer
}

plot_simulations <- function(summary, outdir) {
  if(!requireNamespace("ggplot2",quietly=TRUE)) stop("Install ggplot2 for figures.")
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  plots<-list()
  shapes<-c("PDR-PI1"=16,"PDR-PI2"=17,"PDR-CV1"=15,"PDR-CV2"=18,
    "ProxLL-raw"=1,"ProxLL-US"=0,"ProxLL-raw-matched"=2,"ProxLL-US-matched"=5)
  if(nrow(summary$test)) {
    d<-summary$test
    for(s in intersect(c("CC","CM","MC"),unique(d$scenario))) {
      z<-d[d$scenario==s,]
      p<-ggplot2::ggplot(z,ggplot2::aes(x=.data$delta,y=.data$rejection_probability,
        color=.data$method,shape=.data$method))+ggplot2::geom_line()+ggplot2::geom_point()+
        ggplot2::scale_shape_manual(values=shapes)+ggplot2::facet_grid(null~n)+ggplot2::theme_bw()+
        ggplot2::labs(x="Departure",y="Rejection probability",color="Method",shape="Method")
      plots[[paste0("power_",s)]]<-p
    }
    z<-d[d$delta==0 & d$scenario %in% c("CC","CM","MC"),]
    plots$size<-ggplot2::ggplot(z,ggplot2::aes(x=.data$n,y=.data$rejection_probability,
      color=.data$method,shape=.data$method))+ggplot2::geom_hline(yintercept=.05,color="grey50")+
      ggplot2::geom_point()+ggplot2::scale_shape_manual(values=shapes)+ggplot2::facet_grid(null~scenario)+
      ggplot2::theme_bw()+ggplot2::labs(x="Sample size",y="Rejection probability")
    if(length(unique(z$n))>1L) plots$size<-plots$size+ggplot2::geom_line()
  }
  if(nrow(summary$curve)) {
    for(e in unique(summary$curve$estimand)) {
      z<-summary$curve[summary$curve$estimand==e & summary$curve$scenario %in% c("CC","CM","MC"),]
      plots[[paste0(e,"_coverage")]]<-ggplot2::ggplot(z,
        ggplot2::aes(x=.data$a,y=.data$coverage,color=.data$method))+
        ggplot2::geom_hline(yintercept=.95,color="grey50")+ggplot2::geom_line()+
        ggplot2::facet_grid(scenario~n)+ggplot2::theme_bw()+
        ggplot2::labs(x="Exposure",y="Pointwise coverage")
    }
  }
  for(name in names(plots)) ggplot2::ggsave(file.path(outdir,paste0(name,".png")),
    plots[[name]],width=12,height=8,dpi=180)
  invisible(plots)
}
