#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

static double choose_small(int n, int k) {
  if (k < 0 || k > n) return 0.0;
  double value = 1.0;
  for (int j = 1; j <= k; ++j) value *= (double)(n - k + j) / (double)j;
  return value;
}

static double range_sum(const std::vector<double>& prefix, int n,
                        int power, int left, int right) {
  const int stride = n + 1;
  return prefix[power * stride + right] - prefix[power * stride + left];
}

static double centered_power_sum(const std::vector<double>& prefix, int n,
                                 int left, int right, double center,
                                 double bandwidth, int power) {
  double out = 0.0;
  for (int k = 0; k <= power; ++k) {
    out += choose_small(power, k) * std::pow(-center, power - k) *
      range_sum(prefix, n, k, left, right);
  }
  return out / std::pow(bandwidth, power);
}

static bool cholesky4(const double a[4][4], double l[4][4],
                      double tolerance, double& pivot_ratio) {
  double scale = 1.0;
  for (int i = 0; i < 4; ++i) scale = std::max(scale, std::fabs(a[i][i]));
  for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) l[i][j] = 0.0;
  double min_pivot = R_PosInf, max_pivot = 0.0;
  for (int i = 0; i < 4; ++i) {
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
  return R_finite(pivot_ratio) && pivot_ratio > 0;
}

static void chol_solve4(const double l[4][4], const double b[4], double x[4]) {
  double z[4];
  for (int i = 0; i < 4; ++i) {
    double value = b[i];
    for (int k = 0; k < i; ++k) value -= l[i][k] * z[k];
    z[i] = value / l[i][i];
  }
  for (int i = 3; i >= 0; --i) {
    double value = z[i];
    for (int k = i + 1; k < 4; ++k) value -= l[k][i] * x[k];
    x[i] = value / l[i][i];
  }
}

// Exact O(n) Epanechnikov local-linear and local-quadratic-intercept LOOCV.
// [[Rcpp::export]]
Rcpp::List epan_level_loocv_cpp(Rcpp::NumericVector x,
                                Rcpp::NumericMatrix y,
                                Rcpp::NumericVector bandwidths,
                                double trim_lower,
                                double trim_upper,
                                int min_neighbors = 20,
                                double leverage_limit = 0.95,
                                double determinant_tolerance = 1e-11) {
  const int n = x.size(), s_count = y.ncol(), h_count = bandwidths.size();
  if (y.nrow() != n) stop("x/y dimension mismatch");
  std::vector<int> order(n);
  for (int i = 0; i < n; ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](int i, int j) { return x[i] < x[j]; });
  std::vector<double> xs(n), ys(n * s_count);
  for (int i = 0; i < n; ++i) {
    xs[i] = x[order[i]];
    for (int s = 0; s < s_count; ++s) ys[s * n + i] = y(order[i], s);
  }
  const int stride = n + 1;
  std::vector<double> px(7 * stride, 0.0), py(s_count * 5 * stride, 0.0);
  for (int i = 0; i < n; ++i) {
    double xp = 1.0;
    for (int k = 0; k <= 6; ++k) {
      px[k * stride + i + 1] = px[k * stride + i] + xp; xp *= xs[i];
    }
    for (int s = 0; s < s_count; ++s) {
      xp = 1.0;
      for (int k = 0; k <= 4; ++k) {
        int offset = (s * 5 + k) * stride;
        py[offset + i + 1] = py[offset + i] + ys[s * n + i] * xp;
        xp *= xs[i];
      }
    }
  }
  NumericMatrix risk_ll(h_count, s_count), risk_pdr(h_count, s_count);
  IntegerVector invalid(h_count), used(h_count);
  for (int hi = 0; hi < h_count; ++hi) {
    double h = bandwidths[hi];
    if (!R_finite(h) || h <= 0) {
      for (int s = 0; s < s_count; ++s) risk_ll(hi,s) = risk_pdr(hi,s) = R_PosInf;
      invalid[hi] = n; continue;
    }
    std::vector<double> sum_ll(s_count,0.0), sum_pdr(s_count,0.0);
    std::vector<int> count_ll(s_count,0), count_pdr(s_count,0);
    int left = 0, right = 0, expected = 0;
    for (int i = 0; i < n; ++i) {
      double a = xs[i];
      while (left < n && xs[left] < a - h) ++left;
      if (right < left) right = left;
      while (right < n && xs[right] <= a + h) ++right;
      if (a < trim_lower || a > trim_upper) continue;
      ++expected;
      if (right - left - 1 < min_neighbors) { ++invalid[hi]; continue; }
      double u[7], t[5];
      for (int k = 0; k <= 6; ++k)
        u[k] = centered_power_sum(px,n,left,right,a,h,k);
      for (int k = 0; k <= 4; ++k) t[k] = u[k] - u[k+2];
      double det2 = t[0]*t[2]-t[1]*t[1];
      double scale2 = std::max(std::max(std::fabs(t[0]),std::fabs(t[2])),1.0);
      if (!R_finite(det2) || det2 <= determinant_tolerance*scale2*scale2) {
        ++invalid[hi]; continue;
      }
      double i00=t[2]/det2, i01=-t[1]/det2;
      double aa=t[0],bb=t[1],cc=t[2],dd=t[2],ee=t[3],ff=t[4];
      double det3=aa*(dd*ff-ee*ee)-bb*(bb*ff-cc*ee)+cc*(bb*ee-cc*dd);
      double scale3=std::max(std::max(std::fabs(aa),std::fabs(dd)),
                             std::max(std::fabs(ff),1.0));
      if (!R_finite(det3) || det3 <= determinant_tolerance*scale3*scale3*scale3) {
        ++invalid[hi]; continue;
      }
      double j02=(bb*ee-cc*dd)/det3, j12=(bb*cc-aa*ee)/det3,
        j22=(aa*dd-bb*bb)/det3;
      double c2=i00*t[2]+i01*t[3];
      double lev_ll=i00, lev_pdr=i00-c2*j02;
      if (!R_finite(lev_ll)||!R_finite(lev_pdr)||
          std::fabs(lev_ll)>=leverage_limit||std::fabs(lev_pdr)>=leverage_limit) {
        ++invalid[hi]; continue;
      }
      for (int s=0;s<s_count;++s) {
        double ry[3];
        for (int k=0;k<=2;++k) {
          double centered0=0.0, centered2=0.0;
          for (int q=0;q<=k;++q) {
            int off=(s*5+q)*stride;
            centered0 += choose_small(k,q)*std::pow(-a,k-q)*
              (py[off+right]-py[off+left]);
          }
          centered0 /= std::pow(h,k);
          for (int q=0;q<=k+2;++q) {
            int off=(s*5+q)*stride;
            centered2 += choose_small(k+2,q)*std::pow(-a,k+2-q)*
              (py[off+right]-py[off+left]);
          }
          centered2 /= std::pow(h,k+2);
          ry[k]=centered0-centered2;
        }
        double fit_ll=i00*ry[0]+i01*ry[1];
        double beta2=j02*ry[0]+j12*ry[1]+j22*ry[2];
        double fit_pdr=fit_ll-c2*beta2;
        double res_ll=(ys[s*n+i]-fit_ll)/(1.0-lev_ll);
        double res_pdr=(ys[s*n+i]-fit_pdr)/(1.0-lev_pdr);
        if (R_finite(res_ll)) { sum_ll[s]+=res_ll*res_ll; ++count_ll[s]; }
        if (R_finite(res_pdr)) { sum_pdr[s]+=res_pdr*res_pdr; ++count_pdr[s]; }
      }
    }
    used[hi]=expected;
    for (int s=0;s<s_count;++s) {
      bool enough_ll=count_ll[s]>=std::max(1,(int)std::floor(0.98*expected));
      bool enough_pdr=count_pdr[s]>=std::max(1,(int)std::floor(0.98*expected));
      risk_ll(hi,s)=enough_ll?sum_ll[s]/count_ll[s]:R_PosInf;
      risk_pdr(hi,s)=enough_pdr?sum_pdr[s]/count_pdr[s]:R_PosInf;
    }
  }
  colnames(risk_ll)=colnames(y); colnames(risk_pdr)=colnames(y);
  return List::create(_["bandwidth"]=bandwidths,_["risk_ll"]=risk_ll,
                      _["risk_pdr"]=risk_pdr,_["invalid_evaluations"]=invalid,
                      _["trimmed_evaluations"]=used);
}

// Exact leave-one-out Epanechnikov density at each observed x.
// [[Rcpp::export]]
Rcpp::NumericVector epan_loo_density_cpp(Rcpp::NumericVector x, double bandwidth) {
  const int n=x.size();
  if (n<2 || !R_finite(bandwidth) || bandwidth<=0) stop("Invalid KDE input");
  std::vector<int> order(n);
  for(int i=0;i<n;++i) order[i]=i;
  std::sort(order.begin(),order.end(),[&](int i,int j){return x[i]<x[j];});
  std::vector<double> xs(n);
  for(int i=0;i<n;++i) xs[i]=x[order[i]];
  const int stride=n+1;
  std::vector<double> prefix(3*stride,0.0);
  for(int i=0;i<n;++i){
    double xp=1.0;
    for(int k=0;k<=2;++k){prefix[k*stride+i+1]=prefix[k*stride+i]+xp;xp*=xs[i];}
  }
  NumericVector out(n);
  int left=0,right=0;
  for(int i=0;i<n;++i){
    double a=xs[i];
    while(left<n && xs[left]<a-bandwidth) ++left;
    if(right<left) right=left;
    while(right<n && xs[right]<=a+bandwidth) ++right;
    double u0=centered_power_sum(prefix,n,left,right,a,bandwidth,0)-1.0;
    double u2=centered_power_sum(prefix,n,left,right,a,bandwidth,2);
    double density=0.75*(u0-u2)/(bandwidth*(n-1.0));
    out[order[i]]=density;
  }
  return out;
}

static void taper_value(double a, double outer_lower, double outer_upper,
                        double flat_lower, double flat_upper,
                        double& rho, double& rho_prime) {
  if (a <= outer_lower || a >= outer_upper) { rho=0.0; rho_prime=0.0; return; }
  if (a >= flat_lower && a <= flat_upper) { rho=1.0; rho_prime=0.0; return; }
  double t, direction, width;
  if (a < flat_lower) {
    width=flat_lower-outer_lower; t=(flat_lower-a)/width; direction=-1.0;
  } else {
    width=outer_upper-flat_upper; t=(a-flat_upper)/width; direction=1.0;
  }
  // 1 - 10t^3 + 15t^4 - 6t^5: C2, with value/first derivative zero at t=1.
  rho=1.0-10.0*std::pow(t,3)+15.0*std::pow(t,4)-6.0*std::pow(t,5);
  double derivative_t=-30.0*t*t+60.0*std::pow(t,3)-30.0*std::pow(t,4);
  rho_prime=derivative_t*direction/width;
}

// Tapered least-squares cross-validation for the literal local-linear slope and
// the final b=h debiased derivative operator, which is exactly the degree-three
// local-polynomial slope.  The criterion is an importance-sampling estimate of
// derivative ISE, with the cross term obtained exactly by integration by parts.
// [[Rcpp::export]]
Rcpp::List epan_derivative_loocv_cpp(Rcpp::NumericVector x,
                                     Rcpp::NumericMatrix y,
                                     Rcpp::NumericVector inverse_density,
                                     Rcpp::NumericVector bandwidths,
                                     double outer_lower,
                                     double outer_upper,
                                     double flat_lower,
                                     double flat_upper,
                                     int min_neighbors=20,
                                     double pivot_tolerance=1e-12) {
  const int n=x.size(), s_count=y.ncol(), h_count=bandwidths.size();
  if(y.nrow()!=n || inverse_density.size()!=n) stop("Derivative CV dimension mismatch");
  std::vector<int> order(n);
  for(int i=0;i<n;++i) order[i]=i;
  std::sort(order.begin(),order.end(),[&](int i,int j){return x[i]<x[j];});
  std::vector<double> xs(n), invf(n), ys(n*s_count);
  for(int i=0;i<n;++i){
    xs[i]=x[order[i]]; invf[i]=inverse_density[order[i]];
    for(int s=0;s<s_count;++s) ys[s*n+i]=y(order[i],s);
  }
  NumericMatrix risk_ll(h_count,s_count), risk_pdr(h_count,s_count);
  IntegerVector invalid_ll(h_count), invalid_pdr(h_count), used(h_count),
    min_neighbors_seen(h_count);
  NumericVector min_pivot_ratio(h_count);
  for(int hi=0;hi<h_count;++hi){
    double h=bandwidths[hi];
    if(!R_finite(h)||h<=0){
      for(int s=0;s<s_count;++s)risk_ll(hi,s)=risk_pdr(hi,s)=R_PosInf;
      invalid_ll[hi]=invalid_pdr[hi]=n; min_neighbors_seen[hi]=0;
      min_pivot_ratio[hi]=0; continue;
    }
    std::vector<double> sum_ll(s_count,0.0), sum_pdr(s_count,0.0);
    std::vector<int> count_ll(s_count,0), count_pdr(s_count,0);
    const int first_target=(int)(std::upper_bound(xs.begin(),xs.end(),outer_lower)-xs.begin());
    const int last_target=(int)(std::lower_bound(xs.begin(),xs.end(),outer_upper)-xs.begin());
    const int expected=last_target-first_target;
    if(expected<=0){
      for(int s=0;s<s_count;++s)risk_ll(hi,s)=risk_pdr(hi,s)=R_PosInf;
      invalid_ll[hi]=invalid_pdr[hi]=n;min_neighbors_seen[hi]=0;
      min_pivot_ratio[hi]=0;continue;
    }
    int left=(int)(std::lower_bound(xs.begin(),xs.end(),xs[first_target]-h)-xs.begin());
    int right=(int)(std::upper_bound(xs.begin(),xs.end(),xs[first_target]+h)-xs.begin());
    int min_n=n; double min_ratio=R_PosInf;
    std::vector<long double> running_x(9,0.0L), running_y(s_count*6,0.0L);
    auto rebuild_moments=[&](double center){
      std::fill(running_x.begin(),running_x.end(),0.0L);
      std::fill(running_y.begin(),running_y.end(),0.0L);
      for(int obs=left;obs<right;++obs){
        long double u=((long double)xs[obs]-center)/h, power=1.0L;
        for(int k=0;k<=8;++k){running_x[k]+=power;power*=u;}
        for(int s=0;s<s_count;++s){
          power=1.0L;long double response=ys[s*n+obs];
          for(int k=0;k<=5;++k){running_y[s*6+k]+=response*power;power*=u;}
        }
      }
    };
    auto update_observation=[&](int obs,double center,long double sign){
      long double u=((long double)xs[obs]-center)/h,power=1.0L;
      for(int k=0;k<=8;++k){running_x[k]+=sign*power;power*=u;}
      for(int s=0;s<s_count;++s){
        power=1.0L;long double response=ys[s*n+obs];
        for(int k=0;k<=5;++k){running_y[s*6+k]+=sign*response*power;power*=u;}
      }
    };
    double anchor=xs[first_target];
    const int reanchor_period=std::max(16,std::min(512,n/20));
    rebuild_moments(anchor);
    for(int i=first_target;i<last_target;++i){
      double a=xs[i];
      if(i>first_target){
        while(left<n&&xs[left]<a-h){update_observation(left,anchor,-1.0L);++left;}
        while(right<n&&xs[right]<=a+h){update_observation(right,anchor,1.0L);++right;}
        if((i-first_target)%reanchor_period==0){anchor=a;rebuild_moments(anchor);}
      }
      int neighbors=right-left-1; min_n=std::min(min_n,neighbors);
      if(neighbors<min_neighbors||!R_finite(invf[i])||invf[i]<=0){
        ++invalid_ll[hi];++invalid_pdr[hi];continue;
      }
      double u[9];
      const long double delta=((long double)a-anchor)/h;
      for(int k=0;k<=8;++k){
        long double value=0.0L;
        for(int r=0;r<=k;++r)value+=(long double)choose_small(k,r)*
          std::pow(-delta,k-r)*running_x[r];
        u[k]=(double)value;
      }
      u[0]-=1.0;
      double t[7],tdot[7];
      for(int j=0;j<=6;++j)t[j]=u[j]-u[j+2];
      for(int j=0;j<=6;++j){
        double first=(j==0)?0.0:j*u[j-1];
        tdot[j]=-(first-(j+2)*u[j+1])/h;
      }
      double det2=t[0]*t[2]-t[1]*t[1];
      double scale2=std::max(std::max(std::fabs(t[0]),std::fabs(t[2])),1.0);
      bool valid2=R_finite(det2)&&det2>pivot_tolerance*scale2*scale2;
      double d[4][4],ddot[4][4],l[4][4],ratio;
      for(int row=0;row<4;++row)for(int col=0;col<4;++col){
        d[row][col]=t[row+col];ddot[row][col]=tdot[row+col];
      }
      bool valid4=cholesky4(d,l,pivot_tolerance,ratio);
      if(!valid2)++invalid_ll[hi];
      if(!valid4)++invalid_pdr[hi]; else min_ratio=std::min(min_ratio,ratio);
      if(!valid2&&!valid4)continue;
      double rho,rhop;taper_value(a,outer_lower,outer_upper,flat_lower,flat_upper,rho,rhop);
      for(int s=0;s<s_count;++s){
        double yu[6];
        for(int k=0;k<=5;++k){
          long double value=0.0L;
          for(int r=0;r<=k;++r)value+=(long double)choose_small(k,r)*
            std::pow(-delta,k-r)*running_y[s*6+r];
          yu[k]=(double)value;
        }
        yu[0]-=ys[s*n+i];
        double moment[4],moment_dot[4];
        for(int j=0;j<4;++j){
          moment[j]=yu[j]-yu[j+2];
          double first=(j==0)?0.0:j*yu[j-1];
          moment_dot[j]=-(first-(j+2)*yu[j+1])/h;
        }
        if(valid2){
          double inv00=t[2]/det2, inv01=-t[1]/det2, inv11=t[0]/det2;
          double beta0=inv00*moment[0]+inv01*moment[1];
          double beta1=inv01*moment[0]+inv11*moment[1];
          double rhs0=moment_dot[0]-tdot[0]*beta0-tdot[1]*beta1;
          double rhs1=moment_dot[1]-tdot[1]*beta0-tdot[2]*beta1;
          double beta1_dot=inv01*rhs0+inv11*rhs1;
          double derivative=beta1/h, derivative_dot=beta1_dot/h;
          double contribution=invf[i]*(rho*derivative*derivative+
            2.0*ys[s*n+i]*(rho*derivative_dot+rhop*derivative));
          if(R_finite(contribution)){sum_ll[s]+=contribution;++count_ll[s];}
        }
        if(valid4){
          double beta[4]={0,0,0,0};chol_solve4(l,moment,beta);
          double rhs[4];
          for(int row=0;row<4;++row){rhs[row]=moment_dot[row];
            for(int col=0;col<4;++col)rhs[row]-=ddot[row][col]*beta[col];}
          double beta_dot[4]={0,0,0,0};chol_solve4(l,rhs,beta_dot);
          double derivative=beta[1]/h, derivative_dot=beta_dot[1]/h;
          double contribution=invf[i]*(rho*derivative*derivative+
            2.0*ys[s*n+i]*(rho*derivative_dot+rhop*derivative));
          if(R_finite(contribution)){sum_pdr[s]+=contribution;++count_pdr[s];}
        }
      }
    }
    used[hi]=expected;min_neighbors_seen[hi]=(min_n==n)?0:min_n;
    min_pivot_ratio[hi]=R_finite(min_ratio)?min_ratio:0.0;
    for(int s=0;s<s_count;++s){
      bool enough_ll=count_ll[s]>=std::max(1,(int)std::floor(0.98*expected));
      bool enough_pdr=count_pdr[s]>=std::max(1,(int)std::floor(0.98*expected));
      risk_ll(hi,s)=enough_ll?sum_ll[s]/count_ll[s]:R_PosInf;
      risk_pdr(hi,s)=enough_pdr?sum_pdr[s]/count_pdr[s]:R_PosInf;
    }
  }
  colnames(risk_ll)=colnames(y);colnames(risk_pdr)=colnames(y);
  return List::create(_["bandwidth"]=bandwidths,_["risk_ll"]=risk_ll,
    _["risk_pdr"]=risk_pdr,_["invalid_ll"]=invalid_ll,
    _["invalid_pdr"]=invalid_pdr,_["tapered_evaluations"]=used,
    _["min_neighbors"]=min_neighbors_seen,_["min_pivot_ratio"]=min_pivot_ratio);
}
