// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bw_boot
double bw_boot(NumericVector x, double h, double g);
RcppExport SEXP _baggingbwsel_bw_boot(SEXP xSEXP, SEXP hSEXP, SEXP gSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type g(gSEXP);
    rcpp_result_gen = Rcpp::wrap(bw_boot(x, h, g));
    return rcpp_result_gen;
END_RCPP
}
// Cbw_boot
double Cbw_boot(int n, double d, NumericVector cnt, double h, double g);
RcppExport SEXP _baggingbwsel_Cbw_boot(SEXP nSEXP, SEXP dSEXP, SEXP cntSEXP, SEXP hSEXP, SEXP gSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cnt(cntSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type g(gSEXP);
    rcpp_result_gen = Rcpp::wrap(Cbw_boot(n, d, cnt, h, g));
    return rcpp_result_gen;
END_RCPP
}
// Cbw_ucv_nb
double Cbw_ucv_nb(NumericVector x, double h);
RcppExport SEXP _baggingbwsel_Cbw_ucv_nb(SEXP xSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(Cbw_ucv_nb(x, h));
    return rcpp_result_gen;
END_RCPP
}
// Cbw_ucv
double Cbw_ucv(int n, double d, NumericVector cnt, double h);
RcppExport SEXP _baggingbwsel_Cbw_ucv(SEXP nSEXP, SEXP dSEXP, SEXP cntSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cnt(cntSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(Cbw_ucv(n, d, cnt, h));
    return rcpp_result_gen;
END_RCPP
}
// bw_den
List bw_den(int nbin, NumericVector x);
RcppExport SEXP _baggingbwsel_bw_den(SEXP nbinSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nbin(nbinSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(bw_den(nbin, x));
    return rcpp_result_gen;
END_RCPP
}
// bw_den_binned
NumericVector bw_den_binned(IntegerVector sx);
RcppExport SEXP _baggingbwsel_bw_den_binned(SEXP sxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type sx(sxSEXP);
    rcpp_result_gen = Rcpp::wrap(bw_den_binned(sx));
    return rcpp_result_gen;
END_RCPP
}
// prnw
NumericVector prnw(const NumericVector& x, const NumericVector& y, double h, double x0);
RcppExport SEXP _baggingbwsel_prnw(SEXP xSEXP, SEXP ySEXP, SEXP hSEXP, SEXP x0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    rcpp_result_gen = Rcpp::wrap(prnw(x, y, h, x0));
    return rcpp_result_gen;
END_RCPP
}
// nw_cv
double nw_cv(const NumericVector& x, const NumericVector& y, double h);
RcppExport SEXP _baggingbwsel_nw_cv(SEXP xSEXP, SEXP ySEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(nw_cv(x, y, h));
    return rcpp_result_gen;
END_RCPP
}
// nw_cv_binning
double nw_cv_binning(const NumericVector& x, const NumericVector& y, int nb, double d, double h);
RcppExport SEXP _baggingbwsel_nw_cv_binning(SEXP xSEXP, SEXP ySEXP, SEXP nbSEXP, SEXP dSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(nw_cv_binning(x, y, nb, d, h));
    return rcpp_result_gen;
END_RCPP
}
// new_nw_cv_binning
double new_nw_cv_binning(const NumericVector& x, const NumericVector& y, const IntegerVector& ind, int nb, double d, double h);
RcppExport SEXP _baggingbwsel_new_nw_cv_binning(SEXP xSEXP, SEXP ySEXP, SEXP indSEXP, SEXP nbSEXP, SEXP dSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind(indSEXP);
    Rcpp::traits::input_parameter< int >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(new_nw_cv_binning(x, y, ind, nb, d, h));
    return rcpp_result_gen;
END_RCPP
}
// nw_binning
double nw_binning(int k, const NumericVector& x, const NumericVector& y, int nb, double d, double h);
RcppExport SEXP _baggingbwsel_nw_binning(SEXP kSEXP, SEXP xSEXP, SEXP ySEXP, SEXP nbSEXP, SEXP dSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(nw_binning(k, x, y, nb, d, h));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_baggingbwsel_bw_boot", (DL_FUNC) &_baggingbwsel_bw_boot, 3},
    {"_baggingbwsel_Cbw_boot", (DL_FUNC) &_baggingbwsel_Cbw_boot, 5},
    {"_baggingbwsel_Cbw_ucv_nb", (DL_FUNC) &_baggingbwsel_Cbw_ucv_nb, 2},
    {"_baggingbwsel_Cbw_ucv", (DL_FUNC) &_baggingbwsel_Cbw_ucv, 4},
    {"_baggingbwsel_bw_den", (DL_FUNC) &_baggingbwsel_bw_den, 2},
    {"_baggingbwsel_bw_den_binned", (DL_FUNC) &_baggingbwsel_bw_den_binned, 1},
    {"_baggingbwsel_prnw", (DL_FUNC) &_baggingbwsel_prnw, 4},
    {"_baggingbwsel_nw_cv", (DL_FUNC) &_baggingbwsel_nw_cv, 3},
    {"_baggingbwsel_nw_cv_binning", (DL_FUNC) &_baggingbwsel_nw_cv_binning, 5},
    {"_baggingbwsel_new_nw_cv_binning", (DL_FUNC) &_baggingbwsel_new_nw_cv_binning, 6},
    {"_baggingbwsel_nw_binning", (DL_FUNC) &_baggingbwsel_nw_binning, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_baggingbwsel(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
