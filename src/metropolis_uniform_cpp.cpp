#include <Rcpp.h>
#include "hit_and_run_fun.h"
#include "adaptive_fun.h"
using namespace Rcpp;

// [[Rcpp::export]]
List metropolis_uniform_cpp(
    IntegerVector init, 
    IntegerMatrix moves,
    int iter, 
    int burn,
    int thin, 
    bool hit_and_run, 
    bool adaptive
){
  IntegerVector current = clone(init);
  int n_total_samples = iter * thin;         // total number of steps
  int n_cells= current.size();                  // number of cells
  int n_moves = moves.ncol();               // number of moves
  IntegerMatrix steps(n_cells, iter);            // columns are states
  IntegerVector which_move(n_total_samples);  // move selection
  IntegerVector proposal(n_cells);               // the proposed moves
  double transition_prob;                             // the probability of transition
  int accept_count = 0;
  int total_count = 0;
  bool any_negative_cells;
  IntegerVector move(n_cells);
  double accept_prob = 0;
  
  Function sample("sample");
  which_move = sample(n_moves, burn + n_total_samples, 1);
  
  // enforce burn-in
  if (burn > 0) for (int j = 0; j < burn; ++j) {
    
    // determine move
    for (int k = 0; k < n_cells; ++k) move[k] = moves(k, which_move[j] - 1);
    
    // compute proposal
    if(hit_and_run) proposal = hit_and_run_fun(current, move);
    if(adaptive) proposal = adaptive_fun(current, move);
    if(hit_and_run == false & adaptive == false) {
      for(int k = 0; k < n_cells; ++k) proposal[k] = current[k] + move[k];
    }    
    // compute probability of transition
    any_negative_cells = false;
    for (int k = 0; k < n_cells; ++k) {
      if (proposal[k] < 0) any_negative_cells = true;
    }
    
    // compute transition probability
    if (any_negative_cells) {
      transition_prob = 0.;
    } else {
      transition_prob = 1.;
    }
    
    // make move
    if (runif(1)[0] < transition_prob) {
      for (int k = 0; k < n_cells; ++k) current[k] = proposal[k];
    }
    
  }
  
  // set first step
  for (int k = 0; k < n_cells; ++k) steps(k,0) = current[k];
  
  // main chain
  for(int i = 0; i < iter; ++i){
    
    // cycle through thinning moves, the last of which is kept
    for(int j = 0; j < thin; ++j){
      
      // determine move
      for(int k = 0; k < n_cells; ++k) move[k] = moves(k, which_move[burn+thin*i+j]-1);
      
      // compute proposal
      if(hit_and_run) proposal = hit_and_run_fun(current, move);
      
      if(adaptive) proposal = adaptive_fun(current, move);
      
      if(hit_and_run == false & adaptive == false) {
        for(int k = 0; k < n_cells; ++k) proposal[k] = current[k] + move[k];
      }
      
      // compute probability of transition
      any_negative_cells = false;
      for(int k = 0; k < n_cells; ++k){
        if(proposal[k] < 0) any_negative_cells = true; 
      }
      
      // compute transition probability
      if (any_negative_cells) {
        transition_prob = 0.;
      } else {
        transition_prob = 1.;
      }
      
        // make move
        if(runif(1)[0] < transition_prob){        
          for(int k = 0; k < n_cells; ++k) current[k] = proposal[k];
          accept_count++;
        }
      
      
      // update total_count
      total_count++;
    }
    
    // record current position in steps
    for(int k = 0; k < n_cells; ++k) steps(k,i) = current[k];
    
  }
  
  // create out list
  List out = List::create(
    Rcpp::Named("steps") = steps,
    Rcpp::Named("accept_count") = accept_count,
    Rcpp::Named("total_count") = total_count,
    Rcpp::Named("accept_prob") = (accept_count + 0.) / (total_count + 0.)
  );
  
  return out;
}