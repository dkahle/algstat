#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix walk(IntegerVector current, IntegerMatrix moves, int iter, int thin){

  int nTotalSamples = iter * thin;         // total number of steps
  int n = current.size();                  // number of cells
  int nMoves = moves.ncol();               // number of moves
  IntegerMatrix steps(n, iter);            // columns are states
  IntegerVector whichMove(nTotalSamples);  // move selection
  NumericVector unifs(nTotalSamples);      // for transition probabilities
  IntegerVector proposal(n);               // the proposed moves
  double prob;                             // the probability of transition
  bool anyIsNegative;
  IntegerVector move(n);

  Function sample("sample");
  whichMove = sample(nMoves, nTotalSamples, 1);
  Function runif("runif");
  unifs = runif(nTotalSamples);
  Function print("print");

  for(int i = 0; i < iter; ++i){
    for(int j = 0; j < thin; ++j){

      // make move
      for(int k = 0; k < n; ++k){
        move[k] = moves(k, whichMove[thin*i+j]-1);
      }

      // compute proposal
      for(int k = 0; k < n; ++k){
        proposal[k] = current[k] + move[k];
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
        prob = 1;
      }

      if(prob > 1){
        prob = 1;
      }

      // make move
      if(unifs[thin*i+j] < prob){
        for(int k = 0; k < n; ++k){
          current[k] = proposal[k];
        }
      }

    }

    // assign state move
    for(int k = 0; k < n; ++k){
      steps(k,i) = current[k];
    }
  }

  return steps;
}
