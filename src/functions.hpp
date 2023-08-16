#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::uvec findin(arma::uvec &a, arma::uvec &b);
arma::vec sfun(arma::vec x, double th);
arma::vec scad(arma::vec v, double lam, const double nu = 1 , const double gam = 3 );
arma::mat inverse(const arma::mat &x, const arma::mat &z, const double nu = 1);
arma::mat inversex(const arma::mat &x, const double nu);
arma::mat inverser( arma::vec indexy, arma::mat &x, arma::mat &z, double nu);
arma::mat inverserx( arma::vec indexy, arma::mat &x,double nu);
arma::sp_mat Dfun(int n);
arma::uvec getgroup(const arma::mat &b2pair, const int n, const double tol);

double nsfun(double x, double th);
double nscad(double v, double lam, const double nu =1 ,
             const double gam  = 3);
arma::uvec ngetgroup(arma::vec b2value,
                     int n, double tol);

#endif
