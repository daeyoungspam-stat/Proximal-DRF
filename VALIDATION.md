# Validation record

Validation environment: R 4.6.1, macOS arm64, R-linked numerical libraries.

`R CMD check --no-manual --no-build-vignettes` completed with `Status: OK`,
including task-plan and historical-recovery routing tests.

## Numerical checks

* Analytic curve and derivative identities, bridge representations, polynomial reproduction, equivalent-kernel moments, accelerated cross-validation, influence-function covariance, and deterministic joint execution passed the bundled mathematical validation.
* The local-quadratic derivative extension passed polynomial, bandwidth, and full-pipeline checks.
* Direct bootstrap response generation, null refitting, smoothing, and numerical integration agreed with the matrix implementation for degrees 0, 1, and 2, both residualizations, and all testing methods. The maximum reported discrepancy was approximately 1.14e-14.
* Unequal fold weights, two-point multiplier moments, corrected p-values, and observed-statistic identities passed.
* The BLAS bootstrap path agreed with the matrix reference for statistics, p-values, and rejection decisions.
* Recovery checks passed for sparse-window admissibility, original candidate risks, the raw/US bandwidth relation, preservation of the original PI1 pilot, and original random-number mappings.

## Execution and comparison

* Two test-profile replications of each of the joint, derivative-comparator, and null-testing engines completed: six valid artifacts. Streaming summaries contained two replications per cell. Plot generation completed with all method shapes represented.
* Repeating a completed test plan verified and retained its existing artifacts.
* Production-profile n=2500, replication 1 was recomputed for the joint and derivative-comparator engines. Curves, bands, bandwidths, cross-validation diagnostics, and the compared nuisance diagnostics agreed with the saved final-run artifacts at tolerance 1e-10.
* The production task plan contains 44,000 tasks. The testing route contains 3947 original, 5 first-recovery, 35 second-recovery, and 13 third-recovery identities.

Validation scope consists of these numerical checks, reduced executions, and the specified production-replication comparison. A complete fresh production run remains the computational reproduction step described in README.md.
