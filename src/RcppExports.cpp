// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// SLCC1
Rcpp::List SLCC1(arma::vec indexy, arma::vec& y, arma::mat& z, arma::mat& x, arma::vec& weights, arma::vec& wtilde, arma::mat& betam0, double nu, double gam, double lam, int maxiter, double tolabs, double tolrel);
RcppExport SEXP _SAESLCC_SLCC1(SEXP indexySEXP, SEXP ySEXP, SEXP zSEXP, SEXP xSEXP, SEXP weightsSEXP, SEXP wtildeSEXP, SEXP betam0SEXP, SEXP nuSEXP, SEXP gamSEXP, SEXP lamSEXP, SEXP maxiterSEXP, SEXP tolabsSEXP, SEXP tolrelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type indexy(indexySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type wtilde(wtildeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type betam0(betam0SEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type gam(gamSEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type tolabs(tolabsSEXP);
    Rcpp::traits::input_parameter< double >::type tolrel(tolrelSEXP);
    rcpp_result_gen = Rcpp::wrap(SLCC1(indexy, y, z, x, weights, wtilde, betam0, nu, gam, lam, maxiter, tolabs, tolrel));
    return rcpp_result_gen;
END_RCPP
}
// SLCC2
Rcpp::List SLCC2(arma::vec indexy, arma::vec& y, arma::mat& x, arma::vec& weights, arma::vec& wtilde, arma::mat& betam0, double nu, double gam, double lam, int maxiter, double tolabs, double tolrel);
RcppExport SEXP _SAESLCC_SLCC2(SEXP indexySEXP, SEXP ySEXP, SEXP xSEXP, SEXP weightsSEXP, SEXP wtildeSEXP, SEXP betam0SEXP, SEXP nuSEXP, SEXP gamSEXP, SEXP lamSEXP, SEXP maxiterSEXP, SEXP tolabsSEXP, SEXP tolrelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type indexy(indexySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type wtilde(wtildeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type betam0(betam0SEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type gam(gamSEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type tolabs(tolabsSEXP);
    Rcpp::traits::input_parameter< double >::type tolrel(tolrelSEXP);
    rcpp_result_gen = Rcpp::wrap(SLCC2(indexy, y, x, weights, wtilde, betam0, nu, gam, lam, maxiter, tolabs, tolrel));
    return rcpp_result_gen;
END_RCPP
}
// SLCC3
Rcpp::List SLCC3(arma::vec indexy, arma::vec& y, arma::mat& x, arma::vec& group, arma::vec& weights, arma::vec& wtilde, arma::mat& betam0, double nu, double gam, double lam, int maxiter, double tolabs, double tolrel);
RcppExport SEXP _SAESLCC_SLCC3(SEXP indexySEXP, SEXP ySEXP, SEXP xSEXP, SEXP groupSEXP, SEXP weightsSEXP, SEXP wtildeSEXP, SEXP betam0SEXP, SEXP nuSEXP, SEXP gamSEXP, SEXP lamSEXP, SEXP maxiterSEXP, SEXP tolabsSEXP, SEXP tolrelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type indexy(indexySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type group(groupSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type wtilde(wtildeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type betam0(betam0SEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type gam(gamSEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type tolabs(tolabsSEXP);
    Rcpp::traits::input_parameter< double >::type tolrel(tolrelSEXP);
    rcpp_result_gen = Rcpp::wrap(SLCC3(indexy, y, x, group, weights, wtilde, betam0, nu, gam, lam, maxiter, tolabs, tolrel));
    return rcpp_result_gen;
END_RCPP
}
// cal_initialr
arma::mat cal_initialr(arma::vec indexy, arma::vec& y, arma::mat& z, arma::mat& x, arma::vec& wtilde, double lam0);
RcppExport SEXP _SAESLCC_cal_initialr(SEXP indexySEXP, SEXP ySEXP, SEXP zSEXP, SEXP xSEXP, SEXP wtildeSEXP, SEXP lam0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type indexy(indexySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type wtilde(wtildeSEXP);
    Rcpp::traits::input_parameter< double >::type lam0(lam0SEXP);
    rcpp_result_gen = Rcpp::wrap(cal_initialr(indexy, y, z, x, wtilde, lam0));
    return rcpp_result_gen;
END_RCPP
}
// cal_initialrx
arma::mat cal_initialrx(arma::vec indexy, arma::vec& y, arma::mat& x, arma::vec& wtilde, double lam0);
RcppExport SEXP _SAESLCC_cal_initialrx(SEXP indexySEXP, SEXP ySEXP, SEXP xSEXP, SEXP wtildeSEXP, SEXP lam0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type indexy(indexySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type wtilde(wtildeSEXP);
    Rcpp::traits::input_parameter< double >::type lam0(lam0SEXP);
    rcpp_result_gen = Rcpp::wrap(cal_initialrx(indexy, y, x, wtilde, lam0));
    return rcpp_result_gen;
END_RCPP
}
// getgroup
arma::uvec getgroup(const arma::mat& deltam, const int n, const double tol);
RcppExport SEXP _SAESLCC_getgroup(SEXP deltamSEXP, SEXP nSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type deltam(deltamSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(getgroup(deltam, n, tol));
    return rcpp_result_gen;
END_RCPP
}
// ngetgroup
arma::uvec ngetgroup(arma::vec b2value, int n, double tol);
RcppExport SEXP _SAESLCC_ngetgroup(SEXP b2valueSEXP, SEXP nSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type b2value(b2valueSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(ngetgroup(b2value, n, tol));
    return rcpp_result_gen;
END_RCPP
}
// getorder
arma::uvec getorder(arma::sp_mat& Cmat);
RcppExport SEXP _SAESLCC_getorder(SEXP CmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type Cmat(CmatSEXP);
    rcpp_result_gen = Rcpp::wrap(getorder(Cmat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SAESLCC_SLCC1", (DL_FUNC) &_SAESLCC_SLCC1, 13},
    {"_SAESLCC_SLCC2", (DL_FUNC) &_SAESLCC_SLCC2, 12},
    {"_SAESLCC_SLCC3", (DL_FUNC) &_SAESLCC_SLCC3, 13},
    {"_SAESLCC_cal_initialr", (DL_FUNC) &_SAESLCC_cal_initialr, 6},
    {"_SAESLCC_cal_initialrx", (DL_FUNC) &_SAESLCC_cal_initialrx, 5},
    {"_SAESLCC_getgroup", (DL_FUNC) &_SAESLCC_getgroup, 3},
    {"_SAESLCC_ngetgroup", (DL_FUNC) &_SAESLCC_ngetgroup, 3},
    {"_SAESLCC_getorder", (DL_FUNC) &_SAESLCC_getorder, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_SAESLCC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
