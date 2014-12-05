#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computeNMsCpp(NumericMatrix x, NumericVector exp){

  int ncol = x.ncol();
  int n = x.nrow();
  NumericVector out(ncol);
  double chisq;

  for(int i = 0; i < ncol; ++i){
    chisq = 0;
    for(int j = 0; j < n; ++j){
      if(x(j,i) > 0) chisq += pow(x(j,i) - exp[j], 2) / x(j,i);
    }
    out[i] = chisq;
  }

  return out;
}
