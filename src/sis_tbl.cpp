#include <Rcpp.h>
#include "rcddAPI.h"
using namespace Rcpp;

// [[Rcpp::export]]

IntegerVector sis_tbl(IntegerMatrix A, IntegerVector suff_stats) {
  int w = 0; // counter for while loop
  int d = A.nrow(); // number of sufficient statisitcs
  int r = A.ncol(); // number of cells in the table
  IntegerVector tbl(r); //table to return
  NumericMatrix constr(d+r, r+2); // constraint matrix
  NumericVector objfun(r+1); // objective function
  int p = 0; // counter for error purposes 
  IntegerMatrix work_A(A.nrow(), A.ncol()); // work matrix to edit configuration matrix
  IntegerVector work_suff_stats(suff_stats.size()); // work vector to edit sufficient statistics
  int min, max; // integers to calculate our range
  
  Function print("print");
  
  // Elements needed to run the linear program solver
  LogicalVector first(1);
  first[0] = true;
  LogicalVector second(1);
  second[0] = false;
  CharacterVector solver(1);
  solver = "DualSimplex";
  
  
  // Theoretically, the loop will run until "correct" table is produced
  while(w < 1){

    for(int i = 0; i < A.nrow(); ++i){
      for(int j = 0; j < A.ncol(); ++j){
        work_A(i,j) = A(i,j);
      }
    }
    for(int i = 0; i < suff_stats.size(); ++i) work_suff_stats[i] = suff_stats[i];
   
   
    bool lpsolved = true; //Logical for checking purposes

      for(int i = 0; i<r;++i){
        // creating the objective function
        for(int j = 0; j <  r+1; ++j){
           objfun[j] = (j == i+1) ? -1:0;
        }
        // constructing the constraint matrix
        int z = 2;
        for(int k = 0; k < d+r ; ++k){
          for(int l = 0; l < r + 2; ++l){
           
            if(k < d && l == 0) constr(k,l) = 1; //First column for equalities
            
            if(k < d && l == 1) constr(k,l) = work_suff_stats[k];//Second column for sufficient statistics (rhs)
           
            //Rest of the columns for A matrix
            for(int m = 2; m < r +2; ++m){
              if(k < d && l == m) constr(k,l) = work_A(k,m-2);
             }
            
            if(k >= d && l == 0) constr(k,l) = 0;//First column for inequalities
        
            if(k >= d && l == 1)constr(k,l) = 0;//Second column for right hand side
    
            //Constructing coefficients for each cell to be positive
            for(int m = 2; m < r + 2; ++m){
              if(k >= d && l == m) constr(k,l) = (z == m) ? -1:0;
          }
        }
        if(k >= d) ++z;
      }
      //Running linear program solver to find range of possible values(min, max)
      SEXP out1 = lpcdd_f(constr, objfun, first, solver);
      String solution = VECTOR_ELT(out1, 0);
      //If solution Optimal continue, else break and start over
      if(solution == "Optimal"){
        IntegerVector val = VECTOR_ELT(out1,3);
        min = Rcpp::as<int>(val);
      }else{
        lpsolved = false;
        break;
      }
      SEXP out2 = lpcdd_f(constr, objfun, second, solver);
      String solution2 = VECTOR_ELT(out2, 0);
      //If solution Optimal continue, else break and start over
      if(solution2 == "Optimal"){
        IntegerVector val2 = VECTOR_ELT(out2,3);
        max = Rcpp::as<int>(val2);
      }else{
        lpsolved = false;
        break;
      }
      //Calculate the range and sample from that range to populate the table
          if(min == max - 1) min = max;
          
          IntegerVector range = seq(min,max);
          IntegerVector value = sample(range,1);
          tbl[i] = Rcpp::as<int>(value);
      
      // Update constraints(work_A and work_suff_stats)
      IntegerVector index;
      int y = 0;
      
      //Updating work_A by changing non-zero A elements to zero in the column
      //Keep track of where non-zero elements were
      for(int o = 0; o < d; ++o){
        if(work_A(o,i) != 0){
          work_A(o,i) = 0;
          index[y] = o;
          ++y;
        }
      }
      //Updating work_suff_stats where elements of work_A were changed
      int x = 0;
      for(int p = 0; p < d; ++p){
        if(p == index[x]){
          work_suff_stats[p] = work_suff_stats[p] - tbl[i] * A(p,i);
          ++x;
        }
      }
    }
    // If all linear programs are solved, index w and end the loop
    if(lpsolved == true) ++w;
    
    // If error continues, only let it continue two times
    ++p;
    if(p > 2) break;
  }
  return tbl;
}

