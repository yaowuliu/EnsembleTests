// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Get_chistat
arma::mat Get_chistat(arma::mat X, arma::mat Sigma, int B, int m_select, bool is_Sigma_identity);
RcppExport SEXP _EnsembleTests_Get_chistat(SEXP XSEXP, SEXP SigmaSEXP, SEXP BSEXP, SEXP m_selectSEXP, SEXP is_Sigma_identitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type m_select(m_selectSEXP);
    Rcpp::traits::input_parameter< bool >::type is_Sigma_identity(is_Sigma_identitySEXP);
    rcpp_result_gen = Rcpp::wrap(Get_chistat(X, Sigma, B, m_select, is_Sigma_identity));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EnsembleTests_Get_chistat", (DL_FUNC) &_EnsembleTests_Get_chistat, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_EnsembleTests(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
