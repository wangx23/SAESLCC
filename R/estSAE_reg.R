#' based on common regression models, three type of estimators
#' @param area area id
#' @param y response
#' @param x covariates intercept will be added
#' @param Xbar population mean
#' @param wts sampling weights
#' @param N area size
#' 
#' @return list
#' @export


estSAE_reg <- function(area, y, x, Xbar, wts, N)
{
  ni <- as.numeric(table(area))
  ybari <- as.numeric(by(y, area, mean))
  Nbar <- mean(N)
  wtilde  <- rep(1/Nbar, nivec) * wts ### adjusted weight in the algorithm
  Wd <- diag(wtilde)
  betahat <- solve(t(x) %*% Wd %*% x)%*% t(x) %*% Wd %*% y
  
  muhat1 <- Xbar %*% betahat
  
  ### SR estimator
  est_SR <- muhat1 + as.numeric(by(wts *(y - x %*% betahat), area, sum))/N
  
  ## SYN estimator
  est_syn <- muhat1
  
  # comp
  est_comp <- ni/N * ybari + muhat1 - as.numeric(by(x %*% betahat, area, sum))/N
  
  res <- list(beta = betahat, SR = est_SR, syn = est_syn, comp = est_comp)
  return(res)
}