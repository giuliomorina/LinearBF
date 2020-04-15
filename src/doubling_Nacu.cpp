#include <Rcpp.h>
#include <algorithm>
#include "doubling_envelope.h"
using namespace Rcpp;
using namespace std;

double doubling_alpha_cpp(const int n, const int k, const double eps) {
  return(std::min(2*(double)k/n,1-2*eps));
}

double doubling_beta_cpp(const int n, const int k, const double eps) {
  double C1 = 1/((1-1.0/std::sqrt(2))*eps);
  double C2 = 72/(1-exp(-2*(pow(eps,2))));
  double res = doubling_alpha_cpp(n,k,eps);
  double r1_aux = ((double)k/n-(0.5-3*eps));
  double r2_aux = ((double)k/n-1.0/9);
  if(r1_aux > 0) {
    res += C1*r1_aux*sqrt(2.0/n);
  }
  if(r2_aux > 0) {
    res += C2*r2_aux*exp(-2*pow(eps,2)*n);
  }
  return(res);
}

// [[Rcpp::export]]
int doubling_find_n0(const double eps) {
  return(find_n0(1,eps,2.0,true));
}

// [[Rcpp::export]]
Rcpp::List doubling_Nacu_Peres_cpp(double p, const double eps,
                                   int n0 = -1, const int max_iter = -1,
                                   const double verbose = false) {
  if(eps >= 1.0/8) {
    stop("eps must be smaller than 1/8.");
  }
  return(envelope_alg(p,eps,2.0,1,n0,max_iter,verbose,true));
}