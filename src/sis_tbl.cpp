#include <Rcpp.h>
#include "rcddAPI.h"
using namespace Rcpp;

// [[Rcpp::export]]

IntegerVector sis_tbl(IntegerMatrix A, IntegerVector suff_stats) {
  int w = 0;
  int n = A.nrow();
  int n1 = A.ncol();
  IntegerVector tbl(n1);
  NumericMatrix constr(n+n1, n1+2);
  NumericVector objfun(n1+1);
  IntegerVector maxs(n1);
  while(w < 1){
  IntegerMatrix work = A;
  IntegerVector work2 = suff_stats;
  int min;
  int max;
  LogicalVector first(1);
  first[0] = true;
  LogicalVector second(1);
  second[0] = false;
  CharacterVector solver(1);
  solver = "DualSimplex";
  bool isAnyNegative = false;
  bool lpsolved = true;
    for(int i = 0; i<n1;++i){
      for(int j = 0; j <  n1+1; ++j){
        if(j == i+1){
          objfun[j] = -1;
        }else{
          objfun[j] = 0;
        }
      }
      int z = 2;
      for(int k = 0; k < n+n1 ; ++k){
        for(int l = 0; l < n1 + 2; ++l){
          if(k < n && l == 0){
            constr(k,l) = 1;
          }if(k < n && l == 1){
            constr(k,l) = work2[k];
          }
          for(int m = 2; m < n1 +2; ++m){
            if(k < n && l == m){
              constr(k,l) = work(k,m-2);
            }
          }if(k >= n && l == 0){
            constr(k,l) = 0;
          }if(k >= n && l == 1){
            constr(k,l) = 0;
          }
          for(int m = 2; m < n1 + 2; ++m){
            if(k >= n && l == m){
              if(z == m){
                constr(k,l) = -1;
              }else{
                constr(k,l) = 0;
              }
            }
          }
        }
        if(k >= n){
          ++z;
        }
      }
      
      SEXP out1 = lpcdd_f(constr, objfun, first, solver);
      String solution = VECTOR_ELT(out1, 0);
      if(solution == "Optimal"){
        IntegerVector val = VECTOR_ELT(out1,3);
        min = Rcpp::as<int>(val);
      }else{
        lpsolved = false;
        break;
      }
      SEXP out2 = lpcdd_f(constr, objfun, second, solver);
      String solution2 = VECTOR_ELT(out2, 0);
      if(solution2 == "Optimal"){
        IntegerVector val2 = VECTOR_ELT(out2,3);
        max = Rcpp::as<int>(val2);
      }else{
        lpsolved = false;
        break;
      }
        if(min == max){
          tbl[i] = min;
        } else {
          IntegerVector range = seq(min,max);
          IntegerVector value = sample(range,1);
          tbl[i] = Rcpp::as<int>(value);
        }
      // Update constraints(work and work2)
      IntegerVector index;
      int y = 0;
      
      for(int o = 0; o < n; ++o){
        if(work(o,i) == 1){
          work(o,i) = 0;
          index[y] = o;
          ++y;
        }
      }
      int x = 0;
      for(int p = 0; p < n; ++p){
        if(p == index[x]){
          work2[p] = work2[p] - tbl[i];
          ++x;
        }
      }
    }
    //Check if elements are non-zero
    for(int q = 0;q < n; ++q){
      if(tbl[q] < 0){
        isAnyNegative = true;
      }
    }
    if(isAnyNegative == false && lpsolved == true){
      ++w;
    }
  }
  return tbl;
}



/*** R
tbl <- c(50, 54, 56, 43, 59, 32, 54, 45, 82)
#Configuration matrix of a 3 by 3 matrix
  A <- hmat(c(3,3),1:2)
  suff_stats <- A %*% t(t(tbl))
  sis_tbl(A, suff_stats)
*/
