########## bootstrap ##############
## obj from SLCC_SR fit
# B number of bootstrap 
## three options: "intercept" (model 1), "creg" (model 2), "ccreg" (model 3)
# group is coordinate group for model 3
# output is the B \times m, (m is the number of area) matrix of estimates and B\times m bootstrap mean
# covariates x
bootstrap_SR <-function(obj, area, y, x, Xbar, wts, N, model, 
                        group=NULL, B, lambda, seed = 123456)
{
  
  nivec <- as.numeric(table(area))
  Nbar <- rep(mean(N), length(N))
  wtilde <- rep(1/Nbar, nivec) * wts ### adjusted weight in the algorithm
  m <- length(N)
  uarea <- unique(area)

  ## refitted estimates 
  res_refit <- obj$refit$refit
  eta_est <- res_refit$eta
  beta_est <- as.matrix(res_refit$beta)
  sig2est <- res_refit$sig2

  #generate population
  N0 <- sum(N)
  Ncm <- c(0,cumsum(N))

  YestB <- matrix(0,B, m)
  YbarB <- matrix(0,B, m)
  
  set.seed(seed)
  for(b in 1:B)
  {
    #### get new data y_star
    e_b <- rnorm(N0)*sqrt(sig2est)
    ebar <- aggregate(e_b, by = list(area =rep(1:m,N)), mean)$x
    
    if(model =="intercept")
    {
      Y_bar_star <- Xbar %*% eta_est + beta_est + ebar
      
      subfun <- function(ii)
      {
        eii <- e_b[(Ncm[ii]+1):Ncm[ii+1]]
        indexii <- area == uarea[ii]
        yii <- x[indexii,,drop = FALSE] %*% eta_est +beta_est[ii] +eii[1:nivec[ii]]
        return(yii)
      }
      
      y_star <- unlist(lapply(1:m, subfun))
    }else{
      Y_bar_star <- rowSums(cbind(1,Xbar) * beta_est) +ebar
      subfun <- function(ii)
      {
        eii <- e_b[(Ncm[ii]+1):Ncm[ii+1]]
        indexii <- area == uarea[ii]
        yii <-  cbind(1,x[indexii,,drop = FALSE]) %*% as.matrix(beta_est[ii,]) + eii[1:nivec[ii]]
        return(yii)
      }
      
      y_star <- unlist(lapply(1:m, subfun))
      
    }
    
    # based on new data, get estimate
    resb <- SLCC_SR(area = area, y = y_star, x = x, Xbar = Xbar,
                    wts = wts, N = N, model = model, lambda = lambda, 
                    group = group)
    yestb <- resb$estY
    
    YbarB[b,] <- Y_bar_star
    YestB[b,] <- yestb
    
  }
  out <- list(YestB = YestB, YbarB = YbarB)
  return(out)
}