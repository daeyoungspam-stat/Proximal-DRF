// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <R_ext/BLAS.h>
// R releases before BLAS_INT was exposed use ordinary Fortran INTEGER.
#ifndef BLAS_INT
#define BLAS_INT int
#endif
using namespace Rcpp;

// Ordinary dense multiplication by the BLAS linked to the selected R module.
// No kernel, coefficient, grid, bandwidth, or statistic is approximated.
// This does not replace R's global BLAS or require a particular vendor library.
// [[Rcpp::export]]
NumericMatrix blas_dense_apply_cpp(NumericMatrix smoother,
                                    NumericMatrix response) {
  if (smoother.ncol() != response.nrow())
    stop("Incompatible dense product dimensions");
  NumericMatrix answer(smoother.nrow(), response.ncol());
  const BLAS_INT m = smoother.nrow(), n = response.ncol(), k = smoother.ncol();
  if (m == 0 || n == 0 || k == 0) return answer;
  const char no_transpose = 'N';
  const double alpha = 1.0, beta = 0.0;
  F77_CALL(dgemm)(&no_transpose, &no_transpose, &m, &n, &k, &alpha,
                  smoother.begin(), &m, response.begin(), &k,
                  &beta, answer.begin(), &m FCONE FCONE);
  return answer;
}
