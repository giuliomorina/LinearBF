#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double alpha_cpp(const int n, const int k, const double eps) {
  return(std::min(2*(double)k/n,1-2*eps));
}

// [[Rcpp::export]]
double beta_cpp(const int n, const int k, const double eps) {
  double C1 = 1/((1-1.0/std::sqrt(2))*eps);
  double C2 = 72/(1-exp(-2*(pow(eps,2))));
  double res = alpha_cpp(n,k,eps);
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
double hypergeom_mean(const NumericVector m, const NumericVector n, 
                      const NumericVector Hn, const double eps,
                      bool compute_a) {
  //Compute E[a(m,Hm)|Hn] where m < n and 
  //Hm|Hn is hypergeometric(n,m,Hn)
  // OR compute E[b(m,Hm)|Hn] (depdening on the value of compute_a)
  NumericVector res(1);
  for(int i=0; i<=Hn[0]; i++) {
    if(compute_a) {
      res[0] = res[0] + exp(lchoose(n-m,Hn-i)[0]+lchoose(m,i)[0]-lchoose(n,Hn)[0]+
        log(alpha_cpp(m[0],i,eps)));
    } else {
      res[0] = res[0] + exp(lchoose(n-m,Hn-i)[0]+lchoose(m,i)[0]-lchoose(n,Hn)[0]+
        log(beta_cpp(m[0],i,eps)));
    }
  }
  return(res[0]);
}

// [[Rcpp::export]]
int find_n0(const double eps) {
  //Find the minimum n such that alpha(n,k) and beta(n,k)
  //are always in [0,1]. Notice that alpha(n,k) is always in [0,1].
  int n0 = 1, k;
  bool greater_one = true;
  while(greater_one) {
    greater_one = false;
    for(k=0; k<=n0; k++) {
      if(beta_cpp(n0,k,eps) >1) {
        greater_one = true;
        break;
      }
    }
    n0 = 2*n0; //Increase n0
  }
  return(n0/2);
}

// [[Rcpp::export]]
int doubling_alg(SEXP toss_coin_func, const double eps, 
                 int n0 = -1, const int max_iter = -1) {
  if(eps >= 1.0/8) {
    stop("eps must be smaller than 1/8.");
  }
  if(n0 <= 0) {
    //Find minimum n0 such that alpha and beta are always in [0,1]
    warning("Prespecifying n0 can speed up the code if the function is reused.");
    n0 = find_n0(eps);
  }
  
  //Get the R function to toss a coin
  Rcpp::Function fn_toss_coin = Rcpp::as<Rcpp::Function>(toss_coin_func);
  
  NumericVector G = runif(1,0,1);
  NumericVector L(1, -9.0), U(1, -9.0);
  NumericVector L_tilde(1, -9.0), U_tilde(1, -9.0);
  NumericVector L_star(1, -9.0), U_star(1, -9.0);
  NumericVector diff_tilde(1, -9.0), diff_tilde_prev(1, -9.0);
  NumericVector H(1, -9.0), H_prev(1, -9.0);
  NumericVector aux(1), aux2(1), aux3(1);
  int n_prev, n; //Current and previous number of total tosses
  int iter = 0; //Number of iterations
  
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
    NumericVector tossed_coins = Rcpp::as<NumericVector>(fn_toss_coin(n-n_prev));
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
    L[0] = alpha_cpp(n,H[0],eps);
    U[0] = beta_cpp(n,H[0],eps);
    //Safety check, L and U need to be in [0,1]
    if(L[0] < 0 || U[0] < 0 || L[0] > 1 || U[0] > 1) {
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
      L_star[0] = hypergeom_mean(aux,aux2,aux3,eps,true);
      U_star[0] = hypergeom_mean(aux,aux2,aux3,eps,false);
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
    //Check final result
    if(G[0] <= L_tilde[0]) {
      return(1);
    } else if(G[0] > U_tilde[0]) {
      return(0);
    }
    //Update n and number of iterations
    n_prev = n; n *= 2; iter += 1;
  }
  warning("Number maximum of iterations exceeded.");
  return(NA_INTEGER);
}
