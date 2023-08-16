#include "functions.hpp"
#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;



// [[Rcpp::export]]
arma::mat cal_initialr(arma::vec indexy, arma::vec &y,arma::mat &z, arma::mat &x,
                       arma::vec &wtilde,
                       double lam0 = 0.0001)
{
  int nt = x.n_rows;
  int p = x.n_cols;
  int q = z.n_cols;

  arma::vec uindexy = unique(indexy);
  int n = uindexy.size();

  int n0 = n*p;
  arma::mat Ip =  eye(p,p);

  arma::sp_mat Dmat = Dfun(n);

  //transformation for y,z,x
  arma::uvec indexi;
  arma::mat xm = zeros(nt,p);
  arma::mat zm = zeros(nt,q);
  arma::vec ym = zeros(nt);

  xm = x.each_col() %sqrt(wtilde);
  zm = z.each_col() %sqrt(wtilde);
  ym = y % sqrt(wtilde);
    
  arma::mat ztz = inv(zm.t() * zm) * zm.t();
  arma::mat tempxy =  trans(xm.each_col() % ym  -
    xm.each_col() %(zm*ztz *ym));

  arma::vec Xty(n0);

  for(int i = 0 ; i < n; i++ )
  {
    indexi = find(indexy == uindexy(i));
    Xty(span(i*p,(i+1)*p - 1)) = sum(tempxy.cols(indexi),1);
  }

  arma::mat Xinv = inverser(indexy,xm,zm,lam0);
  arma::mat reg1 = Xinv * Xty;
  arma::mat betam = trans(reshape(reg1, p, n));
  return(betam);
}


// [[Rcpp::export]]
arma::mat cal_initialrx(arma::vec indexy, arma::vec &y, arma::mat &x,
                        arma::vec &wtilde,
                        double lam0 = 0.0001)
{
  int nt = x.n_rows;
  int p = x.n_cols;

  arma::vec uindexy = unique(indexy);
  int n = uindexy.size();

  int n0 = n*p;
  arma::mat Ip =  eye(p,p);

  arma::sp_mat Dmat = Dfun(n);

  //transformation for y,z,x
  arma::uvec indexi;

  arma::mat xm = zeros(nt,p);
  arma::vec ym = zeros(nt);

  
  xm = x.each_col() %sqrt(wtilde);
  ym = y % sqrt(wtilde);

  arma::mat tempxy =  trans(xm.each_col() % ym);

  arma::vec Xty(n0);

  for(int i = 0 ; i < n; i++ )
  {
    indexi = find(indexy == uindexy(i));
    Xty(span(i*p,(i+1)*p - 1)) = sum(tempxy.cols(indexi),1);
  }

  arma::mat Xinv = inverserx(indexy,xm,lam0);
  arma::mat reg1 = Xinv * Xty;
  arma::mat betam = trans(reshape(reg1, p, n));
  return(betam);
}


