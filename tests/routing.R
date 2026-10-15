library(pdrReproduce)
lines<-readLines(system.file("runner.R",package="pdrReproduce"))
start<-which(startsWith(lines,"engine<-if(task$component"))
end<-which(startsWith(lines,"root<-file.path(request$assets"))
end<-end[end>start][1L]-1L
code<-parse(text=lines[start:end])
p<-simulation_plan("parametric_null")
routing<-vapply(seq_len(nrow(p)),function(i) {
  e<-new.env();e$task<-as.list(p[i,]);eval(code,e);e$engine
},character(1))
stopifnot(sum(routing=="null_original")==3947L,
          sum(routing=="null_recovery5")==5L,
          sum(routing=="null_recovery35")==35L,
          sum(routing=="null_recovery13")==13L,
          all(p$n[routing=="null_recovery5"]==5000L),
          all(p$n[routing %in% c("null_recovery35","null_recovery13")]==2500L))
