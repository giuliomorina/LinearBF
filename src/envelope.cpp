#include <Rcpp.h>
#include <algorithm>
#include "doubling_envelope.h"
using namespace Rcpp;
using namespace std;

double alpha_beta_cpp(const int type, bool compute_a, 
                      const int n, const int k, const double eps,
                      const double C) {
  switch(type) {
  case 1:
    // Doubling as in Nacu-Peres
    if((int)C != 2) {
      stop("The algorithm is valid only for C=2");
    }
    return(compute_a ? doubling_alpha_cpp(n,k,eps) : doubling_beta_cpp(n,k,eps));
    break;
  case 2:
    // Doubling as in Flegal.Herbei
    return(compute_a ? doubling_alpha_Flegal_cpp(n,k,eps,C) : doubling_beta_Flegal_cpp(n,k,eps,C));
    break;
  case 3:
    // Doubling as in Flegal-Herbei with newly proposed envelopes
    return(compute_a ? doubling_alpha_Flegal_Morina_cpp(n,k,eps,C) : doubling_beta_Flegal_Morina_cpp(n,k,eps,C));
    break;
  }
  stop("Type not valid");
}

double hypergeom_mean(const NumericVector m, const NumericVector n, 
                      const NumericVector Hn, const double eps,
                      const double C,
                      const int type, const bool compute_a) {
  //Compute E[a(m,Hm)|Hn] where m < n and 
  //Hm|Hn is hypergeometric(n,m,Hn)
  // OR compute E[b(m,Hm)|Hn] (depdening on the value of compute_a)
  NumericVector res(1);
  for(int i=0; i<=Hn[0]; i++) {
    res[0] = res[0] + exp(lchoose(n-m,Hn-i)[0]+lchoose(m,i)[0]-lchoose(n,Hn)[0]+
      log(alpha_beta_cpp(type,compute_a,m[0],i,eps, C)));
  }
  return(res[0]);
}

int find_n0(const int type, const double eps, const double C, const bool double_n) {
  //Find the minimum n such that alpha(n,k) and beta(n,k)
  //are always in [0,1]. 
  
  int n0 = 1, k;
  double alpha, beta;
  bool not_found = true;
  while(not_found) {
    not_found = false;
    for(k=0; k<=n0; k++) {
      alpha = alpha_beta_cpp(type,true,n0,k,eps,C);
      beta = alpha_beta_cpp(type,false,n0,k,eps,C);
      if(alpha < 0 || alpha > 1 || beta < 0 || beta > 1) {
        not_found = true;
        break;
      }
    }
    n0 = double_n ? 2*n0 : n0+1; //Increase n0
  }
  return(double_n ? n0/2 : n0-1);
}

Rcpp::List envelope_alg(double p, const double eps, const double C, const int type,
                 int n0 = -1, const int max_iter = -1,
                 const double verbose = false,
                 const bool double_n = true) {
  if(n0 <= 0) {
    //Find minimum n0 such that alpha and beta are always in [0,1]
    warning("Prespecifying n0 can speed up the code if the function is reused.");
    n0 = find_n0(type, eps, C, double_n);
    if(verbose) Rcout << "Value of n0:" << n0 << "\n";
  }
  
  NumericVector probs = NumericVector::create(p,1.0-p);
  NumericVector G = runif(1,0,1);
  NumericVector L(1, -9.0), U(1, -9.0);
  NumericVector L_tilde(1, -9.0), U_tilde(1, -9.0);
  NumericVector L_star(1, -9.0), U_star(1, -9.0);
  NumericVector diff_tilde(1, -9.0), diff_tilde_prev(1, -9.0);
  NumericVector H(1, -9.0), H_prev(1, -9.0);
  NumericVector aux(1), aux2(1), aux3(1);
  int n_prev, n, num_tosses = 0; //Current, previous and total number of total tosses
  int iter = 0; //Number of iterations
  int res = NA_INTEGER; //Final result
  
  //First iteration: n=0
  L[0] = 0; U[0] = 1;
  L_tilde[0] = 0; U_tilde[0] = 1;
  diff_tilde[0] = 1;  diff_tilde_prev[0] = 1;
  n_prev = 0;
  //All the others elements are not assigned (-9.0)
  
  //From second iteration
  n = n0; //Set initial number of tosses
  iter = 1;
  while(true && (max_iter <= 0 || iter <= max_iter)) {
    //Toss coin
    NumericVector tossed_coins = Rcpp::as<NumericVector>(Rcpp::sample(2, n-n_prev, true, probs));
    num_tosses += n-n_prev;
    if(iter == 1) {
      H[0] = 0;
    } else {
      H[0] = H_prev[0];
    }
    //Count heads (which are equal to 1)
    for(int j=0; j<tossed_coins.size(); j++) {
      if(tossed_coins[j] == 1) {
        H[0]++;
      }
    }
    
    //Compute L and U
    L[0] = alpha_beta_cpp(type,true,n,H[0],eps,C);
    U[0] = alpha_beta_cpp(type,false,n,H[0],eps,C);
    //Safety check, L and U need to be in [0,1]
    if(L[0] < 0 || U[0] < 0 || L[0] > 1 || U[0] > 1 || U[0] < L[0]) {
      Rcout << n0 << " - " << H << " --> " << L << " - " << U << endl;
      stop("Invalid value for L and U. Increase the initial number of tosses.");
    }
    //Compute L_star and U_star
    //Hypergeom mean uses vectors as input
    if(iter == 1) {
      L_star[0] = 0; // E[L0|F1] = E[0|F1] = 0
      U_star[0] = 1;
    } else {
      aux[0] = n_prev; aux2[0] = n, aux3[0] = H[0];
      L_star[0] = hypergeom_mean(aux,aux2,aux3,eps,C,type,true);
      U_star[0] = hypergeom_mean(aux,aux2,aux3,eps,C,type,false);
    }
    //Compute L_tilde and U_tilde
    if(abs(U_star[0] - L_star[0]) > 1e-16) {
      L_tilde[0] = L_tilde[0] + 
        exp(log(L[0]-L_star[0])-log(U_star[0]-L_star[0])+log(diff_tilde_prev[0]));
      U_tilde[0] = U_tilde[0] -
        exp(log(U_star[0]-U[0])-log(U_star[0]-L_star[0])+log(diff_tilde_prev[0]));
      diff_tilde[0] = exp(log(diff_tilde_prev[0])+log(U[0]-L[0])-log(U_star[0]-L_star[0]));
    }
    //Update results
    H_prev[0] = H[0];
    diff_tilde_prev[0] = diff_tilde[0];
    // Verbose
    if(verbose) Rcout << "Iteration #" << iter << ", U = " << G <<
      ", Interval: [" << L[0] << ", " << U[0] << "]" <<
      ", Interval star: [" << L_star[0] << ", " << U_star[0] << "]" <<
      ", Interval tilde: [" << L_tilde[0] << ", " << U_tilde[0] << "] \n";
    //Check final result
    if(G[0] <= L_tilde[0]) {
      res = 1;
      break;
    } else if(G[0] > U_tilde[0]) {
      res = 0;
      break;
    }
    //Update n and number of iterations
    n_prev = n; 
    n = double_n ? 2*n : n+1; 
    iter += 1;
  }
  if(max_iter > 0 && iter > max_iter) warning("Number maximum of iterations exceeded.");
  return Rcpp::List::create(Rcpp::Named("res") = res,
                            Rcpp::Named("num_tosses") = num_tosses);
}
