#include <Rcpp.h>
#include "sis_tbl.h"
using namespace Rcpp;

// [[Rcpp::export]]
List metropolis_uniform_cpp(
    IntegerVector current, 
    IntegerMatrix moves,
    IntegerVector suff_stats,
    IntegerMatrix config,
    int iter, int thin, 
    bool hit_and_run, 
    bool SIS, bool non_uniform
){

  int nTotalSamples = iter * thin;         // total number of steps
  int n = current.size();                  // number of cells
  int nMoves = moves.ncol();               // number of moves
  IntegerMatrix steps(n, iter);            // columns are states
  IntegerVector whichMove(nTotalSamples);  // move selection
  NumericVector unifs(nTotalSamples);      // for transition probabilities
  NumericVector unifs2(nTotalSamples);
  NumericVector unifs3(nTotalSamples);
  IntegerVector proposal(n);               // the proposed moves
  double prob;                             // the probability of transition
  bool anyIsNegative;
  IntegerVector move(n);
  double acceptProb = 0;
  IntegerVector current_num;
  IntegerVector move_num;
  IntegerVector stepSize;
  IntegerVector upperBound;
  IntegerVector lowerBound;
  int lb;
  int ub;
  IntegerVector run;
  
  whichMove = sample(nMoves, nTotalSamples, 1);
  Function runif("runif");
  unifs = runif(nTotalSamples);
  unifs2 = runif(nTotalSamples);
  unifs3 = runif(nTotalSamples);
  Function print("print");
  
  NumericVector move_dist = rep(1.0, nMoves);
  double counter = moves.ncol();
  int which_move;
  
  for(int i = 0; i < iter; ++i){
    for(int j = 0; j < thin; ++j){
      
      if(non_uniform == true){
        for(int l = 0; l < nMoves; ++l){
          double sums = 0;
          for(int m = 0; m < l+1; ++m){
            sums = sums + move_dist[m];
          }
          
          if(unifs3[thin*i+j] <= sums / counter){
            
            for(int k = 0; k < n; ++k){
              move[k] = moves(k, l);
            }
            which_move = l;
            break;
          }
        }
        for(int k = 0; k < n; ++k){
          proposal[k] = current[k] + move[k];
        }
      }else{
        
        // make move
        for(int k = 0; k < n; ++k){
          move[k] = moves(k, whichMove[thin*i+j]-1);
        }
        if(hit_and_run == true){
          current_num = current[move != 0];
          move_num = move[move != 0];
          stepSize = (-1 * current_num) / move_num;
          lowerBound = stepSize[stepSize < 0];
          upperBound = stepSize[stepSize > 0];
          lb = max(lowerBound);
          ub = min(upperBound);
          IntegerVector test1 = current + lb * move;
          IntegerVector test2 = current + ub * move;
          for(int i = 0; i < n; ++i){
            if(test1[i] < 0){
              lb = 1;
            }
            if(test2[i] < 0){
              ub = -1;
            }
          }
          IntegerVector range = seq(lb,ub);
          run = sample(range,1);
          if(run[1] == 0){
            run[1] = 1;
          }
        }
        if(hit_and_run == TRUE){
          for(int k = 0; k < n; ++k){
            proposal[k] = current[k] + as<int>(run) * move[k];
          }
        }else{
          for(int k = 0; k < n; ++k){
            proposal[k] = current[k] + move[k];
          }
        }
      }
      if(SIS){
        if(unifs2[i] < .05){
          proposal = sis_tbl(config, suff_stats);
        }
      }
      // compute probability of transition
      anyIsNegative = false;
      for(int k = 0; k < n; ++k){
        if(proposal[k] < 0){
          anyIsNegative = true;
        }
      }

      if(anyIsNegative){
        prob = 0;
      } else {
        prob = 1; // accept every proposal = uniform
      }

      if(prob > 1){
        prob = 1;
      }

      // store acceptance probability
      acceptProb = acceptProb + prob / nTotalSamples;

      if(non_uniform == true){
        
        if(unifs[thin*i+j] < prob){
          for(int k = 0; k < n; ++k){
            current[k] = proposal[k];
          }
          
          move_dist[which_move] = move_dist[which_move] + 1;
          ++counter;
        }
      }else{
        // make move
        if(unifs[thin*i+j] < prob){
          
          for(int k = 0; k < n; ++k){
            current[k] = proposal[k];
          }
          
        }
      }

    }

    // assign state move
    for(int k = 0; k < n; ++k){
      steps(k,i) = current[k];
    }
  }

  // create out list
  List out = List::create(
    Rcpp::Named("steps") = steps,
    Rcpp::Named("acceptProb") = acceptProb
  );

  return out;
}
