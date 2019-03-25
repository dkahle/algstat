#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List metropolis_uniform_cpp(
    IntegerVector init, 
    IntegerMatrix moves, 
    int iter, 
    int thin
){
  
  IntegerVector current = clone(init);
  int n_total_samples = iter * thin;         // total number of steps
  int n_cells = current.size();              // number of cells
  int n_moves = moves.ncol();                // number of moves
  IntegerMatrix steps(n_cells, iter);        // columns are states
  IntegerVector which_move(n_total_samples); // move selection
  NumericVector unifs(n_total_samples);      // for transition probabilities
  IntegerVector proposal(n_cells);           // the proposed moves
  double prob;                               // the probability of transition
  bool any_negative_cells;
  IntegerVector move(n_cells);
  double accept_prob = 0;
  
  Function sample("sample");
  which_move = sample(n_moves, n_total_samples, 1);
  Function runif("runif");
  unifs = runif(n_total_samples);
  
  
  // set first step
  for (int k = 0; k < n_cells; ++k) steps(k,0) = current[k];
  
  
  // run main chain 
  for (int i = 1; i < iter; ++i) {
    
    // cycle through thinning moves, the last of which is kept
    for (int j = 0; j < thin; ++j) {
      
      // determine move
      for (int k = 0; k < n_cells; ++k) move[k] = moves(k, which_move[thin*i+j]-1);
      
      // compute proposal
      for (int k = 0; k < n_cells; ++k) proposal[k] = current[k] + move[k];
      
      // compute probability of transition
      any_negative_cells = false;
      for (int k = 0; k < n_cells; ++k) {
        if (proposal[k] < 0) any_negative_cells = true;
      }
      
      // compute transition probability
      if (any_negative_cells) {
        prob = 0;
      } else {
        prob = 1;
      }
      
      
      
      
      // store acceptance probability
      accept_prob = accept_prob + prob / n_total_samples;
      
      // make move
      if (unifs[thin*i+j] < prob) {
        for (int k = 0; k < n_cells; ++k) current[k] = proposal[k];
      }
      
    }
    
    // record current position in steps
    for (int k = 0; k < n_cells; ++k) steps(k,i) = current[k];
    
  }
  
  // create out list
  List out = List::create(
    Rcpp::Named("steps") = steps,
    Rcpp::Named("accept_prob") = accept_prob
  );
  
  return out;
}
