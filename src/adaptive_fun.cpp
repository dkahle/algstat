#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector adaptive_fun(IntegerVector current, IntegerVector move) {

  int n = current.size(); 
  IntegerVector current_num;
  IntegerVector move_num;
  IntegerVector stepSize;
  IntegerVector upperBound;
  IntegerVector lowerBound;
  int lb;
  int ub;
  IntegerVector run;
  IntegerVector constant = IntegerVector::create(-1,1);
  int constant2;
  IntegerVector w_current(n);
  IntegerVector w_proposal(n);
  IntegerVector proposal(n);
  bool anyIsNegative;
  double prob; 
  bool didMove;
  NumericVector unifs(1000); 
  Function runif("runif");
  unifs = runif(1000);
  
  current_num = current[move != 0];
  move_num = move[move != 0];
  stepSize = (-1 * current_num) / move_num;
  lowerBound = stepSize[stepSize < 0];
  upperBound = stepSize[stepSize > 0];
  lb = max(lowerBound);
  ub = min(upperBound);
  
  if(is_true(any(stepSize == 0))){
    IntegerVector test1 = current + lb * move;
    IntegerVector test2 = current + ub * move;
    for(int i = 0; i < n; ++i){
      if(test1[i] < 0) lb = 1;
      if(test2[i] < 0) ub = -1;
    }
  }
  
  int line_length = ub-lb + 1;
  if(line_length < 0) line_length = 1;
  
  for(int m = 0; m < n;++m){
    w_current[m] = current[m];
  }
  
  for(int l = 0; l < line_length; ++l){
    
    constant2 = as<int>(Rcpp::sample(constant, 1));
    if(constant2 == 0) constant2 = 1;
    
    for(int k = 0; k < n;++k){
      w_proposal[k] = w_current[k] + constant2 * move[k];
    }
    
    anyIsNegative = false;
    for(int k = 0; k < n; ++k){
      if(w_proposal[k] < 0){
        anyIsNegative = true;
      }
    }
    if(anyIsNegative){
      prob = 0;
    } else {
      prob = exp( sum(lgamma(w_current+1)) - sum(lgamma(w_proposal+1)) );
    }
    
    if(prob > 1){
      prob = 1;
    }
    // make move
    if(unifs[l] < prob) {
      for(int k = 0; k < n; ++k){
        w_current[k] = w_proposal[k];
      }
    }
  }
  
  didMove = false;
  for(int k = 0; k < n; ++k) {
    if(w_current[k] != current[k]) didMove = true;
  }
  if(didMove == true) {
    for(int k = 0; k < n; ++k) {
      proposal[k] = w_current[k];
    }
  } else {
    constant2 = as<int>(Rcpp::sample(constant, 1));
    for(int k = 0; k < n; ++k){
      proposal[k] = current[k] + constant2 * move[k];
    }
  }
  
  return proposal;
}



/*** R
adaptive_fun(c(8,3,8,1), c(1,-1,-1,1))
*/
