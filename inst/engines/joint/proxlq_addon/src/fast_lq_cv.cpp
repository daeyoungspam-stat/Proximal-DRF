#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

static double choose_small_lq(int n, int k) {
  if (k < 0 || k > n) return 0.0;
  double value = 1.0;
  for (int j = 1; j <= k; ++j)
    value *= static_cast<double>(n - k + j) / static_cast<double>(j);
  return value;
}

static bool cholesky3_lq(const double a[3][3], double l[3][3],
                         double tolerance, double& pivot_ratio) {
  double scale = 1.0;
  for (int i = 0; i < 3; ++i) scale = std::max(scale, std::fabs(a[i][i]));
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) l[i][j] = 0.0;
  double min_pivot = R_PosInf, max_pivot = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j <= i; ++j) {
      double value = a[i][j];
      for (int k = 0; k < j; ++k) value -= l[i][k] * l[j][k];
      if (i == j) {
        if (!R_finite(value) || value <= tolerance * scale) return false;
        l[i][j] = std::sqrt(value);
        min_pivot = std::min(min_pivot, value);
        max_pivot = std::max(max_pivot, value);
      } else {
        l[i][j] = value / l[j][j];
      }
    }
  }
  pivot_ratio = min_pivot / max_pivot;
  return R_finite(pivot_ratio) && pivot_ratio > 0.0;
}

static void chol_solve3_lq(const double l[3][3], const double b[3], double x[3]) {
  double z[3];
  for (int i = 0; i < 3; ++i) {
    double value = b[i];
    for (int k = 0; k < i; ++k) value -= l[i][k] * z[k];
    z[i] = value / l[i][i];
  }
  for (int i = 2; i >= 0; --i) {
    double value = z[i];
    for (int k = i + 1; k < 3; ++k) value -= l[k][i] * x[k];
    x[i] = value / l[i][i];
  }
}

static void taper_value_lq(double a, double outer_lower, double outer_upper,
                           double flat_lower, double flat_upper,
                           double& rho, double& rho_prime) {
  if (a <= outer_lower || a >= outer_upper) {
    rho = 0.0; rho_prime = 0.0; return;
  }
  if (a >= flat_lower && a <= flat_upper) {
    rho = 1.0; rho_prime = 0.0; return;
  }
  double t, direction, width;
  if (a < flat_lower) {
    width = flat_lower - outer_lower;
    t = (flat_lower - a) / width;
    direction = -1.0;
  } else {
    width = outer_upper - flat_upper;
    t = (a - flat_upper) / width;
    direction = 1.0;
  }
  rho = 1.0 - 10.0 * std::pow(t, 3) + 15.0 * std::pow(t, 4) -
    6.0 * std::pow(t, 5);
  const double derivative_t = -30.0 * t * t + 60.0 * std::pow(t, 3) -
    30.0 * std::pow(t, 4);
  rho_prime = derivative_t * direction / width;
}

// Exact O(n) tapered derivative-risk CV for the degree-two local-polynomial
// slope. The local fit and its derivative with respect to the evaluation
// coordinate are both analytic. Epanechnikov's multiplicative constant cancels
// from the local-polynomial coefficient.
// [[Rcpp::export]]
Rcpp::List epan_derivative_lq_loocv_cpp(Rcpp::NumericVector x,
                                        Rcpp::NumericMatrix y,
                                        Rcpp::NumericVector inverse_density,
                                        Rcpp::NumericVector bandwidths,
                                        double outer_lower,
                                        double outer_upper,
                                        double flat_lower,
                                        double flat_upper,
                                        int min_neighbors = 20,
                                        double pivot_tolerance = 1e-12) {
  const int n = x.size(), s_count = y.ncol(), h_count = bandwidths.size();
  if (y.nrow() != n || inverse_density.size() != n)
    stop("Degree-two derivative CV dimension mismatch");

  std::vector<int> order(n);
  for (int i = 0; i < n; ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](int i, int j) { return x[i] < x[j]; });
  std::vector<double> xs(n), invf(n), ys(n * s_count);
  for (int i = 0; i < n; ++i) {
    xs[i] = x[order[i]];
    invf[i] = inverse_density[order[i]];
    for (int s = 0; s < s_count; ++s) ys[s * n + i] = y(order[i], s);
  }

  NumericMatrix risk_lq(h_count, s_count);
  IntegerVector invalid_lq(h_count), used(h_count), min_neighbors_seen(h_count);
  NumericVector min_pivot_ratio(h_count);

  for (int hi = 0; hi < h_count; ++hi) {
    const double h = bandwidths[hi];
    if (!R_finite(h) || h <= 0.0) {
      for (int s = 0; s < s_count; ++s) risk_lq(hi, s) = R_PosInf;
      invalid_lq[hi] = n;
      min_neighbors_seen[hi] = 0;
      min_pivot_ratio[hi] = 0.0;
      continue;
    }

    std::vector<double> risk_sum(s_count, 0.0);
    std::vector<int> risk_count(s_count, 0);
    const int first_target = static_cast<int>(
      std::upper_bound(xs.begin(), xs.end(), outer_lower) - xs.begin());
    const int last_target = static_cast<int>(
      std::lower_bound(xs.begin(), xs.end(), outer_upper) - xs.begin());
    const int expected = last_target - first_target;
    if (expected <= 0) {
      for (int s = 0; s < s_count; ++s) risk_lq(hi, s) = R_PosInf;
      invalid_lq[hi] = n;
      min_neighbors_seen[hi] = 0;
      min_pivot_ratio[hi] = 0.0;
      continue;
    }

    int left = static_cast<int>(
      std::lower_bound(xs.begin(), xs.end(), xs[first_target] - h) - xs.begin());
    int right = static_cast<int>(
      std::upper_bound(xs.begin(), xs.end(), xs[first_target] + h) - xs.begin());
    int minimum_neighbors = n;
    double minimum_ratio = R_PosInf;
    std::vector<long double> running_x(7, 0.0L);
    std::vector<long double> running_y(s_count * 5, 0.0L);

    auto rebuild_moments = [&](double center) {
      std::fill(running_x.begin(), running_x.end(), 0.0L);
      std::fill(running_y.begin(), running_y.end(), 0.0L);
      for (int obs = left; obs < right; ++obs) {
        long double u = (static_cast<long double>(xs[obs]) - center) / h;
        long double power = 1.0L;
        for (int k = 0; k <= 6; ++k) {
          running_x[k] += power;
          power *= u;
        }
        for (int s = 0; s < s_count; ++s) {
          power = 1.0L;
          const long double response = ys[s * n + obs];
          for (int k = 0; k <= 4; ++k) {
            running_y[s * 5 + k] += response * power;
            power *= u;
          }
        }
      }
    };

    auto update_observation = [&](int obs, double center, long double sign) {
      const long double u = (static_cast<long double>(xs[obs]) - center) / h;
      long double power = 1.0L;
      for (int k = 0; k <= 6; ++k) {
        running_x[k] += sign * power;
        power *= u;
      }
      for (int s = 0; s < s_count; ++s) {
        power = 1.0L;
        const long double response = ys[s * n + obs];
        for (int k = 0; k <= 4; ++k) {
          running_y[s * 5 + k] += sign * response * power;
          power *= u;
        }
      }
    };

    double anchor = xs[first_target];
    const int reanchor_period = std::max(16, std::min(512, n / 20));
    rebuild_moments(anchor);

    for (int i = first_target; i < last_target; ++i) {
      const double a = xs[i];
      if (i > first_target) {
        while (left < n && xs[left] < a - h) {
          update_observation(left, anchor, -1.0L); ++left;
        }
        while (right < n && xs[right] <= a + h) {
          update_observation(right, anchor, 1.0L); ++right;
        }
        if ((i - first_target) % reanchor_period == 0) {
          anchor = a;
          rebuild_moments(anchor);
        }
      }

      const int neighbors = right - left - 1;
      minimum_neighbors = std::min(minimum_neighbors, neighbors);
      if (neighbors < min_neighbors || !R_finite(invf[i]) || invf[i] <= 0.0) {
        ++invalid_lq[hi];
        continue;
      }

      double u[7];
      const long double delta = (static_cast<long double>(a) - anchor) / h;
      for (int k = 0; k <= 6; ++k) {
        long double value = 0.0L;
        for (int r = 0; r <= k; ++r)
          value += static_cast<long double>(choose_small_lq(k, r)) *
            std::pow(-delta, k - r) * running_x[r];
        u[k] = static_cast<double>(value);
      }
      u[0] -= 1.0;

      double t[5], t_dot[5];
      for (int j = 0; j <= 4; ++j) t[j] = u[j] - u[j + 2];
      for (int j = 0; j <= 4; ++j) {
        const double first = (j == 0) ? 0.0 : j * u[j - 1];
        t_dot[j] = -(first - (j + 2) * u[j + 1]) / h;
      }

      double d[3][3], d_dot[3][3], l[3][3], ratio = 0.0;
      for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
          d[row][col] = t[row + col];
          d_dot[row][col] = t_dot[row + col];
        }
      }
      if (!cholesky3_lq(d, l, pivot_tolerance, ratio)) {
        ++invalid_lq[hi];
        continue;
      }
      minimum_ratio = std::min(minimum_ratio, ratio);

      double rho, rho_prime;
      taper_value_lq(a, outer_lower, outer_upper, flat_lower, flat_upper,
                     rho, rho_prime);

      for (int s = 0; s < s_count; ++s) {
        double yu[5];
        for (int k = 0; k <= 4; ++k) {
          long double value = 0.0L;
          for (int r = 0; r <= k; ++r)
            value += static_cast<long double>(choose_small_lq(k, r)) *
              std::pow(-delta, k - r) * running_y[s * 5 + r];
          yu[k] = static_cast<double>(value);
        }
        yu[0] -= ys[s * n + i];

        double moment[3], moment_dot[3];
        for (int j = 0; j < 3; ++j) {
          moment[j] = yu[j] - yu[j + 2];
          const double first = (j == 0) ? 0.0 : j * yu[j - 1];
          moment_dot[j] = -(first - (j + 2) * yu[j + 1]) / h;
        }

        double beta[3] = {0.0, 0.0, 0.0};
        chol_solve3_lq(l, moment, beta);
        double rhs[3];
        for (int row = 0; row < 3; ++row) {
          rhs[row] = moment_dot[row];
          for (int col = 0; col < 3; ++col)
            rhs[row] -= d_dot[row][col] * beta[col];
        }
        double beta_dot[3] = {0.0, 0.0, 0.0};
        chol_solve3_lq(l, rhs, beta_dot);

        const double derivative = beta[1] / h;
        const double derivative_dot = beta_dot[1] / h;
        const double contribution = invf[i] *
          (rho * derivative * derivative + 2.0 * ys[s * n + i] *
           (rho * derivative_dot + rho_prime * derivative));
        if (R_finite(contribution)) {
          risk_sum[s] += contribution;
          ++risk_count[s];
        }
      }
    }

    used[hi] = expected;
    min_neighbors_seen[hi] = (minimum_neighbors == n) ? 0 : minimum_neighbors;
    min_pivot_ratio[hi] = R_finite(minimum_ratio) ? minimum_ratio : 0.0;
    for (int s = 0; s < s_count; ++s) {
      const bool enough = risk_count[s] >=
        std::max(1, static_cast<int>(std::floor(0.98 * expected)));
      risk_lq(hi, s) = enough ? risk_sum[s] / risk_count[s] : R_PosInf;
    }
  }

  colnames(risk_lq) = colnames(y);
  return List::create(
    _["bandwidth"] = bandwidths,
    _["risk_lq"] = risk_lq,
    _["invalid_lq"] = invalid_lq,
    _["tapered_evaluations"] = used,
    _["min_neighbors"] = min_neighbors_seen,
    _["min_pivot_ratio"] = min_pivot_ratio
  );
}
