#include <Rcpp.h>
#include <algorithm>
#include "doubling_envelope.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double modified_double_func_cpp(double p, double C, double eps) {
  // This is a modified target function for a linear BF which is twice
  // differentiable, as proposed by Flegal, Herbei (2012) in "Exact sampling
  // for intractable probability distributions via a Bernoulli Factory".
  // The function assume that Cp < 1-eps and fix delta such that 0<delta<omega
  // It is then defined as:
  // a*p                       if p \in [0,(1-eps)/a)
  // (1-eps)+F(p-(1-eps)/a)    if p  in [(1-eps)/a,1]
  // where F(p) is given by:
  // F(p) = delta*\int_0^{a*p/delta} e^{-t^2} dt
  if(p >= 0 && p<=(1.0-eps)/C) {
    return(C*p);
  } else if(p <= 1) {
    return((1-eps) + eps*((R::pnorm(C*(p-(1-eps)/C)/eps,0.0,sqrt(0.5),true,false) - 0.5)*sqrt(M_PI)));
  } else {
    return(NA_REAL);
  }
  
}

double doubling_alpha_Flegal_cpp(const int n, const int k, 
                                 const double eps, const double C) {
  return(modified_double_func_cpp((double)k/n, C, eps));
}

double doubling_beta_Flegal_cpp(const int n, const int k, 
                                 const double eps, const double C) {
  double M = C*C*sqrt(2.0)/(eps*sqrt(M_E));
  double res = doubling_alpha_Flegal_cpp(n,k,eps,C);
  res += M/(2.0*n);
  return(res);
}

// [[Rcpp::export]]
double const_a_n_Flegal_Morina_cpp(const int n, const int k, const double eps,
                                   const double C) {
  // Compute constant bounding second derivative
  double M = C*C*sqrt(2.0)/(eps*sqrt(M_E));
  // Compute the coefficients a_n
  double norm_const = sinh(M_PI*sqrt(M))/(sqrt(M)*M_PI);
  double log_a_n = 0;
  for(int j=2; j<=n; j++) {
    log_a_n += log(1+M/pow((j-1),2));
  }
  log_a_n -= log(norm_const);
  return(exp(log_a_n));
}


double doubling_alpha_Flegal_Morina_cpp(const int n, const int k, const double eps,
                                       const double C) {
  // Compute the bound
  double a_n = const_a_n_Flegal_Morina_cpp(n,k,eps,C);
  double res = a_n*modified_double_func_cpp((double)k/n, C, eps);
  return(res);
}

double doubling_beta_Flegal_Morina_cpp(const int n, const int k, const double eps,
                                       const double C) {
  // Compute the bound
  double a_n = const_a_n_Flegal_Morina_cpp(n,k,eps,C);
  double res = 1-a_n*modified_double_func_cpp((double)k/n, C, eps);
  return(res);
}

// [[Rcpp::export]]
int doubling_Flegal_find_n0(const double eps, const double C) {
  return(find_n0(2,eps,C,false));
}

// [[Rcpp::export]]
Rcpp::List doubling_Flegal_Herbei_cpp(double p, const double eps, const double C,
                                   int n0 = -1, const int max_iter = -1,
                                   const double verbose = false) {
  return(envelope_alg(p,eps,C,2,n0,max_iter,verbose,false));
}


// [[Rcpp::export]]
int doubling_Flegal_Morina_find_n0(const double eps, const double C) {
  return(find_n0(3,eps,C,false));
}

// [[Rcpp::export]]
Rcpp::List doubling_Flegal_Herbei_Morina_cpp(double p, const double eps, const double C,
                                      int n0 = -1, const int max_iter = -1,
                                      const double verbose = false) {
  return(envelope_alg(p,eps,C,3,n0,max_iter,verbose,false));
}

