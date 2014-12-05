#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computeUProbsCpp(NumericMatrix x){

  int ncol = x.ncol();
  int n = x.nrow();
  NumericVector out(ncol);
  double ph = 0;

  for(int i = 0; i < ncol; ++i){
    ph = 0;
    for(int k = 0; k < n; ++k){
      ph = ph + lgamma(x(k,i)+1);
    }
    out[i] = -1 * ph; //exp(0 - ph);
  }

  return out;
}
