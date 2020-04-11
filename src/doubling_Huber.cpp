#include <Rcpp.h>
#include <stack> 
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List logistic_cpp(float C, float p) {
  // f(p) = Cp/(1+Cp)
  // The function is taken from "Designing perfect simulation algorithms
  // using local correctness" by M. Huber (2019), in an iterative form rather than
  // recursive. This is equivalent to the 2-coin algorithm proposed in
  // "Barker's algorithm for Bayesian inference with intractable likelihoods"
  // by F. Goncalves et al. (2017)
  int B, X, res;
  int num_tosses = 0;
  NumericVector probs_logistic = NumericVector::create(C/(1.0+C),1.0-C/(1.0+C));
  NumericVector probs = NumericVector::create(p,1.0-p);

  while(true) {
    B = (int)Rcpp::sample(2, 1, true, probs_logistic)[0];
    if(B != 1) {
      res = 0; break;
    } else {
      X = (int)Rcpp::sample(2, 1, true, probs)[0];
      num_tosses += 1;
      if(X == 1) {
        res = 1; break;
      }
    }
  }
  if(res != 1) res = 0;
  return Rcpp::List::create(Rcpp::Named("res") = res,
                            Rcpp::Named("num_tosses") = num_tosses);
}

// [[Rcpp::export]]
Rcpp::List doubling_Huber_2019_iter_cpp(float C, int i, float eps, float p) {
  // https://www.codeproject.com/Articles/418776/How-to-replace-recursive-functions-using-stack-and
  // f(p) = (Cp)^i, for Cp < 1-eps
  // The function is taken from "Designing perfect simulation algorithms
  // using local correctness" by M. Huber (2019).
  // This is a translation in an iterative form with a custom stack.
  struct DoublingStackStruct {
    float C;        
    int i;
    float eps;
    int stage;  
  };
  // Initialise variables
  float beta;
  NumericVector probs_beta(2);
  int B1, B2;
  Rcpp::List aux;
  // Prepare return variable
  int retVal = -1;
  int num_tosses = 0;
  // Initialise stack
  std::stack<DoublingStackStruct> DoublingStack;
  DoublingStackStruct currentSnapshot;
  currentSnapshot.C = C;
  currentSnapshot.i = i;
  currentSnapshot.eps = eps;
  currentSnapshot.stage = 0; // Initial stage
  DoublingStack.push(currentSnapshot);
  // Start loop
  while(!DoublingStack.empty()) {
    currentSnapshot = DoublingStack.top();
    DoublingStack.pop();
    switch(currentSnapshot.stage) {
    case 0: 
      // Main algorithm
      if(currentSnapshot.i == 0) {
        retVal = 1;
        continue;
      } else if(currentSnapshot.i > 3.55/currentSnapshot.eps) {
        beta = (1.0-currentSnapshot.eps/2)/(1.0-currentSnapshot.eps);
        probs_beta[0] = pow(beta,-currentSnapshot.i);
        probs_beta[1] = 1.0-probs_beta[0];
        B1 = (int)Rcpp::sample(2, 1, true, probs_beta)[0];
        if(B1 != 1) {
          retVal = 0;
          continue;
        } else {
          // Push to the stack
          currentSnapshot.stage = 1;
          DoublingStack.push(currentSnapshot); //Re-add to stack
          // Create new snapshot to add to stack
          DoublingStackStruct newSnapshot;
          newSnapshot.C = beta*currentSnapshot.C;
          newSnapshot.i = currentSnapshot.i;
          newSnapshot.eps = currentSnapshot.eps/2.0;
          newSnapshot.stage = 0;
          DoublingStack.push(newSnapshot);
          continue;
        }
      } else {
        aux = logistic_cpp(currentSnapshot.C, p);
        num_tosses += (int)aux["num_tosses"];
        B2 = (int)aux["res"];
        if(B2 != 1) B2=0;
        // Push to the stack
        currentSnapshot.stage = 2;
        DoublingStack.push(currentSnapshot); //Re-add to stack
        // Create new snapshot to add to stack
        DoublingStackStruct newSnapshot;
        newSnapshot.C = currentSnapshot.C;
        newSnapshot.i = currentSnapshot.i + 1 - 2*B2;
        newSnapshot.eps = currentSnapshot.eps;
        newSnapshot.stage = 0;
        DoublingStack.push(newSnapshot);
        continue;
      }
      break;
    case 1:
      //This contains things that should be done after the 1st recursive call. In
      //this case there are none.
      break;
    case 2:
      //This contains things that should be done after the 2nd recursive call. In
      //this case there are none.
      break;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("res") = retVal,
                            Rcpp::Named("num_tosses") = num_tosses);
}
