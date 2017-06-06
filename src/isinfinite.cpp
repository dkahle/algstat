#include <Rcpp.h>
#include "isinfinite.h"
using namespace Rcpp;


// [[Rcpp::export]]
LogicalVector isinfinite(NumericVector x){
  return is_infinite(x);
}


