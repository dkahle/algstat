#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rfiberOne(IntegerMatrix A, IntegerVector b){

  int m = A.nrow();         // number of suff stats
  int n = A.ncol();         // length of b, x (table)
  IntegerVector x(n);       // table
  IntegerVector bounds(2);  // number of cells
  int L, U;                 // lower and upper bounds
  double prob = 1.0;        // the probability of transition
  bool reject = false;      //
  
  
  // make a copy of b called b2; b2 will be changed
  // without changing the value of b in R
  IntegerVector b2(b.size());
  for(int i = 0; i < b.size(); ++i) b2[i] = b[i];
  //IntegerVector b2 = clone(b); // a bit slower
  

  Function sample("sample.int");
  Function cell_bounds("cell_bounds");
  
  
  // do main columns
  for(int i = 0; i < n; ++i){
  
    // compute cell bounds
    bounds = cell_bounds(A, b2, "lp", i+1); 
    L = bounds[0]; U = bounds[1]; 
    
    // sample cell; note NA_integer_ gets parsed to INT_MIN in cpp
    if(L == INT_MIN || U == INT_MIN || (L > U)){
      x[i]   = -1;
      prob   = 0;
      reject = true;
    } else if(L == U){
      x[i] = L;
    } else {
      x[i] = as<IntegerVector>(sample(U-(L-1), 1, 1))[0] + (L-1);
      prob = prob * 1/(U - L + 1); // = prob * 1/length(L:U)
    }
    
    if(i == n-1 || reject == true) break;
    
    // update cell totals
    for(int r = 0; r < m; ++r){      
      b2[r] = b2[r] - x[i]*A(r,i);
    }        
    
  }
  
  
  // short-circuit if rejecting
  if(reject == true){
    return List::create(
      _["x"] = x,
      _["prob"] = prob,
      _["reject"] = reject
    );
  }
  
  
  // check sampled table
  int sumChecker = 0;
  if(!reject){
    for(int i = 0; i < m; ++i){
      sumChecker = 0;
      for(int j = 0; j < n; ++j){
        sumChecker += A(i,j)*x[j];
      }
      if(sumChecker != b[i]){
        prob   = 0;
        reject = true;
        break;
      }
    }
  }
  
  
  // create out list
  return List::create(
    _["x"] = x,
    _["prob"] = prob,
    _["reject"] = reject
  );
}

