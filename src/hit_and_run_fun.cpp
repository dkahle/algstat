#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]


IntegerVector hit_and_run_fun(IntegerVector current, IntegerVector move) {
  
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
  IntegerVector proposal(n);
  
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
  if(lb > ub){
    run = Rcpp::sample(constant, 1);
    
  } else {
    IntegerVector range = seq(lb,ub);
    run = Rcpp::sample(range,1);
  }
  if(run[0] == 0){
    run = Rcpp::sample(constant, 1);
  }
  
  for(int k = 0; k < n; ++k){
    proposal[k] = current[k] + as<int>(run) * move[k];
  }
  return proposal;
}



