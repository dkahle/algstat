// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// computeCRsCpp
NumericVector computeCRsCpp(NumericMatrix x, NumericVector exp, double lambda);
RcppExport SEXP _algstat_computeCRsCpp(SEXP xSEXP, SEXP expSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type exp(expSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(computeCRsCpp(x, exp, lambda));
    return rcpp_result_gen;
END_RCPP
}
// computeG2sCpp
NumericVector computeG2sCpp(NumericMatrix x, NumericVector exp);
RcppExport SEXP _algstat_computeG2sCpp(SEXP xSEXP, SEXP expSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type exp(expSEXP);
    rcpp_result_gen = Rcpp::wrap(computeG2sCpp(x, exp));
    return rcpp_result_gen;
END_RCPP
}
// computeNMsCpp
NumericVector computeNMsCpp(NumericMatrix x, NumericVector exp);
RcppExport SEXP _algstat_computeNMsCpp(SEXP xSEXP, SEXP expSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type exp(expSEXP);
    rcpp_result_gen = Rcpp::wrap(computeNMsCpp(x, exp));
    return rcpp_result_gen;
END_RCPP
}
// computeUProbsCpp
NumericVector computeUProbsCpp(NumericMatrix x);
RcppExport SEXP _algstat_computeUProbsCpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(computeUProbsCpp(x));
    return rcpp_result_gen;
END_RCPP
}
// computeX2sCpp
NumericVector computeX2sCpp(NumericMatrix x, NumericVector exp);
RcppExport SEXP _algstat_computeX2sCpp(SEXP xSEXP, SEXP expSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type exp(expSEXP);
    rcpp_result_gen = Rcpp::wrap(computeX2sCpp(x, exp));
    return rcpp_result_gen;
END_RCPP
}
// metropolis_hypergeometric_cpp
List metropolis_hypergeometric_cpp(IntegerVector init, IntegerMatrix moves, int iter, int thin);
RcppExport SEXP _algstat_metropolis_hypergeometric_cpp(SEXP initSEXP, SEXP movesSEXP, SEXP iterSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type init(initSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type moves(movesSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(metropolis_hypergeometric_cpp(init, moves, iter, thin));
    return rcpp_result_gen;
END_RCPP
}
// metropolis_uniform_cpp
List metropolis_uniform_cpp(IntegerVector init, IntegerMatrix moves, int iter, int thin);
RcppExport SEXP _algstat_metropolis_uniform_cpp(SEXP initSEXP, SEXP movesSEXP, SEXP iterSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type init(initSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type moves(movesSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(metropolis_uniform_cpp(init, moves, iter, thin));
    return rcpp_result_gen;
END_RCPP
}
// rfiberOne
List rfiberOne(IntegerMatrix A, IntegerVector b);
RcppExport SEXP _algstat_rfiberOne(SEXP ASEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(rfiberOne(A, b));
    return rcpp_result_gen;
END_RCPP
}
// walk
IntegerMatrix walk(IntegerVector current, IntegerMatrix moves, int iter, int thin);
RcppExport SEXP _algstat_walk(SEXP currentSEXP, SEXP movesSEXP, SEXP iterSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type current(currentSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type moves(movesSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(walk(current, moves, iter, thin));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_algstat_computeCRsCpp", (DL_FUNC) &_algstat_computeCRsCpp, 3},
    {"_algstat_computeG2sCpp", (DL_FUNC) &_algstat_computeG2sCpp, 2},
    {"_algstat_computeNMsCpp", (DL_FUNC) &_algstat_computeNMsCpp, 2},
    {"_algstat_computeUProbsCpp", (DL_FUNC) &_algstat_computeUProbsCpp, 1},
    {"_algstat_computeX2sCpp", (DL_FUNC) &_algstat_computeX2sCpp, 2},
    {"_algstat_metropolis_hypergeometric_cpp", (DL_FUNC) &_algstat_metropolis_hypergeometric_cpp, 4},
    {"_algstat_metropolis_uniform_cpp", (DL_FUNC) &_algstat_metropolis_uniform_cpp, 4},
    {"_algstat_rfiberOne", (DL_FUNC) &_algstat_rfiberOne, 2},
    {"_algstat_walk", (DL_FUNC) &_algstat_walk, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_algstat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
