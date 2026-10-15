library(pdrReproduce)
p<-simulation_plan()
stopifnot(nrow(p)==44000L, !anyDuplicated(p$key),
  sum(p$component=="null")==4000L,
  sum(p$component=="joint")==20000L,
  sum(p$component=="proxlq")==20000L)
stopifnot(inherits(try(simulation_plan(reps=c(1,1)),silent=TRUE),"try-error"))
