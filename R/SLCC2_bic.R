#' this code gives regression coefficients for SLCC2 based on BIC
#' 
#' 
#' @description for a given lambda vector, it gives the bic loop and the estimates of regression coef

#' @param area area id
#' @param y response
#' @param x covariates, need to add intercept if needed
#' @param wts sampling weights
#' @param lambda tuning vector
#' 
#' @return list
#' @export

SLCC2_bic <- function(area, y, x, wts, lambda)
{
  nt <- length(y)
  m <- length(unique(area))
  nivec <- as.numeric(table(area))
  ncx <- ncol(x)
  N <- as.numeric(by(wts, area, sum))
  Nbar <- mean(N)
  wtilde  <- rep(1/Nbar, nivec) * wts ### adjusted weight in the algorithm
  weights <- rep(1, m*(m-1)/2) ### pairwise weights (cij) in the algorithm,
  ## can be adjusted for spatial data

  betam01 <- SAESLCC::cal_initialrx(area,y,x,wtilde = wtilde)
  
  nlam <- length(lambda)
  sig2est <- rep(0, nlam)
  Cm <- log(m*ncx)
  
  ngest = rep(0, nlam)
  
  for(j in 1:nlam)
  {
    resj <- SLCC2(area,y = y, x = x,
                  weights = weights, wtilde = wtilde,betam0 = betam01,
                  lam = lambda[j])
    sig2est[j] <- resj$sig2
    ngest[j] <- length(unique(resj$cluster[,1]))
  }
  BICm <- data.frame(sig2 = sig2est, ng = ngest,
                     bic = log(sig2est) + Cm*ngest*ncx*log(nt)/nt)
  
  ress <- SLCC2(area,y = y, x = x,
                weights = weights, wtilde = wtilde,betam0 = betam01,
                lam = lambda[which.min(BICm$bic)])
  

  
  out <- list(fit = ress, bic = BICm, lams = lambda[which.min(BICm$bic)])
  return(out)
}