#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computeX2sCpp(NumericMatrix x, NumericVector exp){

  int ncol = x.ncol();
  int n = x.nrow();
  NumericVector out(ncol);
  double chisq;

  for(int i = 0; i < ncol; ++i){
    chisq = 0;
    for(int j = 0; j < n; ++j){
      if(exp[j] > 0) chisq += pow(x(j,i) - exp[j], 2) / exp[j];
    }
    out[i] = chisq;
  }

  return out;
}
