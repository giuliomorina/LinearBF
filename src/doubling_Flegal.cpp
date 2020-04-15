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

double const_a_n_Flegal_Morina_cpp(const int n, const int k, const double eps,
                                   const double C, const double a_n_m_1) {
  // Compute constant bounding second derivative
  double M = C*C*sqrt(2.0)/(eps*sqrt(M_E));
  // Compute coefficient a_2
  if(n <= 1) {
    return(-9.0);
  }
  if(n == 2) {
    double norm_const = sinh(M_PI*sqrt(M))/(sqrt(M)*M_PI);
    return((1.0+M)/norm_const);
  }
  if(a_n_m_1 <= 0 || a_n_m_1 > 1) {
    stop("a_n_m_1 cannot be negative or greater than 1.");
  }
  // Compute the coefficient a_n = a_n_m_1 * (1+M/(n-1)^2)
  double log_a_n = log(a_n_m_1) + log(1.0+M/pow((n-1),2));
  return(exp(log_a_n));
}


double doubling_alpha_Flegal_Morina_cpp(const int n, const int k, const double eps,
                                       const double C, const double a_n) {
  // Compute the bound
  if(n <= 1) return(0.0);
  double res = a_n*modified_double_func_cpp((double)k/n, C, eps);
  return(res);
}

double doubling_beta_Flegal_Morina_cpp(const int n, const int k, const double eps,
                                       const double C, const double a_n) {
  // Compute the bound
  if(n <= 1) return(1.0);
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

