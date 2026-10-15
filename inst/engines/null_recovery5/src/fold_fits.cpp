#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// Observation-location fits. Local windows are evaluated directly in centered
// coordinates; no subtractive raw-moment prefix sums and no n by n matrix are
// formed. The only zero-fit convention is exact local rank deficiency, detected
// by the number of distinct locations having strictly positive kernel weight.

static long double symmetric_condition(const long double input[3][3], int p) {
  long double a[3][3] = {};
  for (int i = 0; i < p; ++i)
    for (int j = 0; j < p; ++j) a[i][j] = input[i][j];
  const long double eps = std::numeric_limits<long double>::epsilon();
  for (int iteration = 0; iteration < 80; ++iteration) {
    int row = 0, col = 1;
    long double largest = 0.0L, scale = 0.0L;
    for (int i = 0; i < p; ++i) {
      scale = std::max(scale, std::fabs(a[i][i]));
      for (int j = i + 1; j < p; ++j)
        if (std::fabs(a[i][j]) > largest) {
          largest = std::fabs(a[i][j]); row = i; col = j;
        }
    }
    if (largest <= eps * scale) break;
    const long double tau = (a[col][col] - a[row][row]) / (2.0L * a[row][col]);
    const long double t = (tau >= 0.0L ? 1.0L : -1.0L) /
      (std::fabs(tau) + std::sqrt(1.0L + tau * tau));
    const long double c = 1.0L / std::sqrt(1.0L + t * t), s = t * c;
    const long double ar = a[row][row], ac = a[col][col], cross = a[row][col];
    a[row][row] = ar - t * cross;
    a[col][col] = ac + t * cross;
    a[row][col] = a[col][row] = 0.0L;
    for (int k = 0; k < p; ++k) if (k != row && k != col) {
      const long double kr = a[k][row], kc = a[k][col];
      a[k][row] = a[row][k] = c * kr - s * kc;
      a[k][col] = a[col][k] = s * kr + c * kc;
    }
  }
  long double lo = a[0][0], hi = a[0][0];
  for (int i = 1; i < p; ++i) {
    lo = std::min(lo, a[i][i]); hi = std::max(hi, a[i][i]);
  }
  if (!(lo > 0.0L)) return std::numeric_limits<long double>::infinity();
  return hi / lo;
}

static bool solve_positive_system(const long double gram[3][3],
                                  const long double rhs[3], int p,
                                  long double beta[3]) {
  long double l[3][3] = {}, z[3] = {};
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j <= i; ++j) {
      long double value = gram[i][j];
      for (int k = 0; k < j; ++k) value -= l[i][k] * l[j][k];
      if (i == j) {
        if (!(value > 0.0L)) return false;
        l[i][j] = std::sqrt(value);
      } else l[i][j] = value / l[j][j];
    }
    long double value = rhs[i];
    for (int k = 0; k < i; ++k) value -= l[i][k] * z[k];
    z[i] = value / l[i][i];
  }
  for (int i = p - 1; i >= 0; --i) {
    long double value = z[i];
    for (int k = i + 1; k < p; ++k) value -= l[k][i] * beta[k];
    beta[i] = value / l[i][i];
  }
  return std::isfinite(beta[0]);
}

// [[Rcpp::export]]
Rcpp::List epan_observation_fits_cpp(Rcpp::NumericVector a,
                                    Rcpp::NumericVector y,
                                    double h, int degree,
                                    double condition_limit) {
  const int n = a.size(), p = degree + 1;
  if (n != y.size() || n < 1 || degree < 1 || degree > 2 ||
      !std::isfinite(h) || h <= 0.0 || !std::isfinite(condition_limit) ||
      condition_limit <= 1.0) Rcpp::stop("Invalid local-fit inputs");
  std::vector<int> order(n);
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(a[i]) || !std::isfinite(y[i]))
      Rcpp::stop("Nonfinite local-fit inputs");
    order[i] = i;
  }
  std::sort(order.begin(), order.end(), [&](int i, int j) { return a[i] < a[j]; });
  Rcpp::NumericVector fitted(n), condition(n), solve_condition(n);
  Rcpp::IntegerVector local_n(n), unique_local_n(n);
  std::vector<int> zero_index, bad_index;
  int left = 0, right = 0;
  for (int r = 0; r < n; ++r) {
    if ((r & 127) == 0) Rcpp::checkUserInterrupt();
    const int original = order[r];
    const long double center = a[original], bandwidth = h;
    while (left < n && (static_cast<long double>(a[order[left]]) - center) / bandwidth <= -1.0L) ++left;
    if (right < left) right = left;
    while (right < n && (static_cast<long double>(a[order[right]]) - center) / bandwidth < 1.0L) ++right;
    long double gram[3][3] = {}, rhs[3] = {};
    int distinct = 0;
    for (int j = left; j < right; ++j) {
      const int index = order[j];
      const long double u = (static_cast<long double>(a[index]) - center) / bandwidth;
      const long double weight = 1.0L - u * u;
      if (!(weight > 0.0L)) continue;
      ++local_n[original];
      if (j == left || a[index] != a[order[j - 1]]) ++distinct;
      const long double powers[3] = {1.0L, u, u * u};
      for (int i = 0; i < p; ++i) {
        rhs[i] += weight * powers[i] * y[index];
        for (int k = 0; k < p; ++k)
          gram[i][k] += weight * powers[i] * powers[k];
      }
    }
    unique_local_n[original] = distinct;
    if (distinct < p) {
      fitted[original] = 0.0;
      condition[original] = R_PosInf;
      solve_condition[original] = R_PosInf;
      zero_index.push_back(original + 1);
      continue;
    }
    const long double kappa = symmetric_condition(gram, p);
    condition[original] = static_cast<double>(kappa);
    // A diagonal change of polynomial coordinates does not change the fitted
    // curve. It removes conditioning caused solely by all local u being small.
    // The constant column is unchanged, so the intercept is still beta[0].
    long double scaled[3][3] = {}, scaled_rhs[3] = {}, scale[3] = {};
    for (int i = 0; i < p; ++i) scale[i] = std::sqrt(gram[i][i] / gram[0][0]);
    for (int i = 0; i < p; ++i) {
      scaled_rhs[i] = rhs[i] / scale[i];
      for (int j = 0; j < p; ++j) scaled[i][j] = gram[i][j] / (scale[i] * scale[j]);
    }
    const long double scaled_kappa = symmetric_condition(scaled, p);
    solve_condition[original] = static_cast<double>(scaled_kappa);
    long double beta[3] = {};
    if (!std::isfinite(scaled_kappa) || scaled_kappa > condition_limit ||
        !solve_positive_system(scaled, scaled_rhs, p, beta)) {
      fitted[original] = NA_REAL;
      bad_index.push_back(original + 1);
    } else fitted[original] = static_cast<double>(beta[0]);
  }
  std::sort(zero_index.begin(), zero_index.end());
  std::sort(bad_index.begin(), bad_index.end());
  return Rcpp::List::create(Rcpp::Named("fitted") = fitted,
                            Rcpp::Named("condition") = condition,
                            Rcpp::Named("solve_condition") = solve_condition,
                            Rcpp::Named("local_n") = local_n,
                            Rcpp::Named("unique_local_n") = unique_local_n,
                            Rcpp::Named("zero_index") = zero_index,
                            Rcpp::Named("bad_index") = bad_index);
}
