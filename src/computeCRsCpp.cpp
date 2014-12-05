#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computeCRsCpp(NumericMatrix x, NumericVector exp, double lambda){

  int ncol = x.ncol();
  int n = x.nrow();
  NumericVector out(ncol);
  double cr;

  for(int i = 0; i < ncol; ++i){
    cr = 0;
    for(int j = 0; j < n; ++j){
      if(x(j,i) > 0){
        if(exp[j] > 0) cr += x(j,i) * (pow(x(j,i)/exp[j], lambda) - 1); // 0 ow, contributes nothing to sum
      }
    }
    out[i] = (2/(lambda*(lambda+1))) * cr;
  }

  return out;
}
