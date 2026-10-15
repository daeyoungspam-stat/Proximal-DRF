\name{pdrReproduce}
\alias{simulation_plan}
\alias{run_simulations}
\alias{simulation_status}
\alias{summarize_simulations}
\alias{validate_engines}
\alias{plot_simulations}
\title{Reproduce proximal dose-response simulation experiments}
\description{Task planning, frozen-engine validation, resumable execution, streaming summaries, and figures for curve and derivative inference and parametric-null tests.}
\usage{
simulation_plan(study = c("all", "inference", "parametric_null"),
                profile = c("production", "test"), n = NULL, reps = NULL)
run_simulations(outdir, study = c("all", "inference", "parametric_null"),
                profile = c("production", "test"), workers = 1L,
                n = NULL, reps = NULL)
simulation_status(outdir)
summarize_simulations(outdir)
validate_engines(outdir = tempfile("pdr-validation-"))
plot_simulations(summary, outdir)
}
\arguments{
\item{study}{Experiment collection. Inference comprises the joint curve/derivative engine and the derivative comparison engine.}
\item{profile}{Production settings or a reduced computational test.}
\item{n}{Optional integer vector of sample sizes.}
\item{reps}{Optional integer vector of replication identifiers, such as \code{1:10}.}
\item{outdir}{Output directory. Resume with the same study, profile, sample sizes and replication identifiers. Worker count may change.}
\item{workers}{Concurrent independent single-thread R processes.}
\item{summary}{Object returned by \code{summarize_simulations}.}
}
\value{A task table, validation list, summary list, or invisible list of plots. Execution writes checksummed RDS artifacts and task logs.}
\details{Production inference uses 5000 replications per sample size. Production testing uses 1000 replications per sample size and 4999 two-point bootstrap draws. The 53 historical recovery identities select their corresponding original-seed recovery engines. Summary tables select the primary residualization and the log10-undersmoothed local-quadratic derivative comparator.}
\examples{
p <- simulation_plan("parametric_null", "test", reps = 1L)
stopifnot(nrow(p) == 1L)
}
