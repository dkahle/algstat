#include <Rcpp.h>
#include <random>
#include <iostream>
#include "isinfinite.h"
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
List metropolis_uniform_cpp(
    IntegerVector current, 
    IntegerMatrix moves, 
    int iter, int thin,
    bool hit_and_run
){

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
  double acceptProb = 0;
  NumericVector stepSize_1(n);
  NumericVector current_num(n);
  NumericVector move_num(n);
  NumericVector s(nMoves);
  for(int i = 0; i < nMoves; ++i){
    int num = 0;
    for(int j = 0; j < n;++j){
      if(moves(i,j) != 0){
        ++num;
      }
    }
    s[i] = num;
  }
  int ss = min(s);
  NumericVector stepSize(ss);
  NumericVector upperBound(ss);
  NumericVector lowerBound(ss);
  double lb;
  double ub;
  int run;

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
      if(hit_and_run == TRUE){
        current_num = as<NumericVector>(current);
        move_num = as<NumericVector>(move);
        for(int i = 0; i < n; ++i){
          stepSize_1[i] = -1 * current_num[i] / move_num[i];
        }
        stepSize = stepSize_1[isinfinite(stepSize_1) == FALSE];
        lowerBound = stepSize[stepSize < 0];
        upperBound = stepSize[stepSize > 0];
        lb = max(lowerBound);
        ub = min(upperBound);
        //run = sample(seq(lb,ub),1);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dis(lb,ub);
        run = dis(gen);
        if(run == 0){
          run = 1;
        }
      }
      if(hit_and_run == TRUE){
        for(int k = 0; k < n; ++k){
          proposal[k] = current[k] + run * move[k];
        }
      }else{
        for(int k = 0; k < n; ++k){
          proposal[k] = current[k] + move[k];
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

  // create out list
  List out = List::create(
    Rcpp::Named("steps") = steps,
    Rcpp::Named("acceptProb") = acceptProb
  );

  return out;
}
