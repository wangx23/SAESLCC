####### get area estimates based on regression models ####
# obj is the object of SLCC
# Xbar population mean of X for each area
# Zbar population mean of Z for each area
# wts sampling weights
# N population area size
estSAE <- function(obj, indexy, y, Xbar, Zbar = NULL, wts, N)
{
  
  uindexy <- unique(indexy)
  yhat <- obj$yhat
  
  subfun <- function(ii)
  {
    indexii <- indexy == uindexy[ii]
    yii <- y[indexii]
    yhatii <- yhat[indexii]
    
    return(sum((yii - yhatii)*wts[indexii]))
  }
  estp1 <- unlist(lapply(1:length(uindexy),subfun))
  
  eta_est <- obj$eta
  beta_est <- obj$beta
  
  if(is.null(Zbar))
  {
    estY <- rowSums(Xbar * beta_est) + estp1/N
  }else{
    estY <- Zbar %*% eta_est + rowSums(Xbar * beta_est) + estp1/N
  }
  return(as.numeric(estY))
}