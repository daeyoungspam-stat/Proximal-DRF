# Proximal dose-response simulation reproduction

An R package for the final curve-inference, derivative-inference, and parametric-null experiments.

## Installation

R 4.2 or later, a C++ compiler, and the BLAS/Fortran runtime matching R support the numerical engines. Linux and macOS support parallel task execution. Windows supports single-worker execution with the R toolchain.

```r
install.packages(c("Rcpp", "Matrix", "ggplot2"))
```

From the cloned repository directory:

```sh
R CMD INSTALL .
```

## Validation and small execution

```r
library(pdrReproduce)
validate_engines("validation-output")
run_simulations("results/smoke", profile = "test", reps = 1L, workers = 2L)
simulation_status("results/smoke")
s <- summarize_simulations("results/smoke")
s$curve
s$test
plot_simulations(s, "results/smoke/figures")
```

The test profile reduces sample size, bootstrap draws, grids, and replication count. Production settings specify the manuscript experiment. A one-replication execution checks the numerical pipeline; Monte Carlo comparisons require the production replication counts.

## Full production reproduction

```r
run_simulations("results/final", profile = "production", workers = 12L)
simulation_status("results/final")
s <- summarize_simulations("results/final")
plot_simulations(s, "results/final/figures")
```

The default plan contains 20,000 joint curve/derivative tasks, 20,000 derivative-comparator tasks, and 4,000 parametric-null tasks. Each task produces one validated artifact. Runtime and storage depend on the numerical libraries and hardware; the full plan is a substantial computational experiment.

Separate experiment directories support independent execution:

```r
run_simulations("results/inference", study = "inference", workers = 12L)
run_simulations("results/testing", study = "parametric_null", workers = 12L)
```

For selected original replication identities:

```r
run_simulations("results/subset", study = "parametric_null",
                n = 2500L, reps = c(818L, 819L), workers = 2L)
```

## Statistical specification

* Curve and derivative inference share the sinusoidal outcome design, five-fold bridge fitting, sample sizes 2500, 5000, 7500, and 10000, and 5000 replications per sample size. PDR-PI1, PDR-PI2, PDR-CV1, and PDR-CV2 accompany ProxLL-US for the curve and ProxLQ-US for the derivative. Bands use 10,000 Gaussian draws and a 61-point evaluation grid.
* Parametric-null testing covers constant, linear, and quadratic nulls, a common cubic departure, and the grid 0, .004, .008, .012, .016, .024, .032, .040, .052, .064, .080. Each sample size has 1000 replications. The fold-averaged statistic uses fold-size weights, two-point multipliers, 4999 draws, and 512 integration nodes. Raw and log10-undersmoothed local-linear comparators include canonical and matched quadratic-null fits.
* Raw artifacts retain both residualization strategies and the supplementary derivative undersmoothing comparison. The summary API selects `fold_parametric` and the log10 ProxLQ-US result, matching the final primary presentation. Scenario MM remains available as a diagnostic.

## Final-run provenance

The joint inference engine corresponds to configuration `4e7d28a1`; the derivative comparison corresponds to `0bebf0fb`. The final testing aggregate combines 3947 original artifacts (`0754d03f`) with 5, 35, and 13 original-seed recovery artifacts (`0e6839ae`, `1775d37c`, `15981712`). The package routes those 53 replication identities to their respective frozen recovery engines. These engines check admissibility of the original final bandwidth candidates across the required fold operators and preserve the original PI1 pilot rule. This mapping reproduces the completed experiment rather than defining a new universal bandwidth rule.

`inst/original_engine_manifest.rds` records the original MD5 for each selected engine file; `inst/engine_manifest.rds` records the packaged MD5. Numerical estimation and bootstrap routines retain their original source. Three joint-engine files receive fingerprint-only portability changes: relative paths use literal substring extraction, and the derivative extension validates the packaged parent identity. This handles installation paths containing spaces and parentheses. The portable package supplies its own orchestration and a reduced set of execution support files, so packaged configuration fingerprints differ from historical deployment fingerprints. Each run records the packaged source manifest and numerical configuration. Floating-point results depend on R, BLAS, compiler, and package versions; session information accompanies each replication.

## Results and recovery

`parts/` contains checksummed result envelopes; `logs/` contains per-task diagnostics; `checkpoints/` contains completed testing-condition prefixes. Repeating the same call validates existing artifacts and executes remaining tasks with their original seeds. Completed curve and derivative replications are restart boundaries. Testing checkpoints additionally preserve completed condition/scenario columns. An interrupted controller can leave `controller.lock`; remove that directory after confirming that its processes have stopped. Preserve artifacts and checkpoints for diagnosis when a task reports an error.

`summarize_simulations()` requires a complete valid task plan. It returns pointwise bias, empirical standard deviation, RMSE, coverage, interval width, band coverage, and testing rejection probabilities. The saved `summary.rds` contains these tables and the task audit. Plotting produces diagnostic PNG figures from the summary tables.

## Repository metadata

Package metadata contains a project designation and a reserved example-domain maintenance address. The address is a packaging placeholder. The repository contains simulation software and configuration; generated datasets and result directories belong in separate storage.

The LICENSE file reserves rights pending an explicit public reuse-license selection. Installation and numerical validation are recorded in VALIDATION.md.
