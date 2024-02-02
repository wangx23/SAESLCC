####### this code gives estimator for SR ####
## it implies several procedures ####
# step 1: use SLCC function to find group based on standardized data (using a loop)
# step 2: use the group result to construct the new estimator
# area: area id
# lambda, a sequence of lambda
# add intercept automatically

### three models, "intercept" (model 1), "creg" (model 2) and "ccreg" (model 3)
SLCC_SR <- function(area, y, x, Xbar, 
                    wts, N, model, group=NULL, lambda)
{
  nt <- length(y)
  m <- length(unique(area))
  nivec <- as.numeric(table(area))
  wtilde  <- rep(1/N, nivec) * wts ### adjusted weight in the algorithm
  weights <- rep(1, m*(m-1)/2) ### pairwise weights (cij) in the algorithm,
  ## can be adjusted for spatial data
  
  ## ### scale version for model fitting
  y_sc <- scale(y) 
  x_sc <- apply(x,2,scale) ### scale version for model fitting
  xintercept <- matrix(1, nrow = nt)
  
  if(model == "intercept")
  {
    betam01 <- cal_initialr(area,y_sc,x_sc,xintercept,wtilde = wtilde)
  }else{
    betam01 <- cal_initialrx(area,y_sc,cbind(1,x_sc),wtilde = wtilde)
  }
  
  nlam <- length(lambda)
  sig2est <- rep(0, nlam)
  Cm <- log(m)
  
  if(model == "intercept")
  {
    ngest = rep(0, nlam)

    for(j in 1:nlam)
    {
      resj <- SLCC1(indexy = area,y = y_sc,z = x_sc,x = xintercept,
                     weights = weights, wtilde = wtilde,betam0 = betam01,lam = lambda[j])
      sig2est[j] <- resj$sig2
      ngest[j] <- length(unique(resj$cluster[,1]))
    }
    BICm <- data.frame(sig2 = sig2est, ng = ngest) %>% 
      mutate(bic = log(sig2) + Cm*ng*log(nt)/nt)
    
    resm <- SLCC1(indexy = area,y = y_sc,z = x_sc,x = xintercept,
                  weights = weights, wtilde = wtilde,betam0 = betam01,
                  lam = lambda[which.min(BICm$bic)])
    
    refit_SR <- estSAE_SR(obj = resm,indexy = area,y = y,
                           x = xintercept, Xbar = matrix(1,nrow=m), 
                           z = x,Zbar = Xbar,wts = wts, N = N, 
                           model ="intercept")
 
  }
  
  if(model == "creg")
  {
    ngest = rep(0, nlam)
    
    for(j in 1:nlam)
    {
      resj <- SLCC2(indexy = area,y = y_sc, x = cbind(1, x_sc),
                    weights = weights, wtilde = wtilde,betam0 = betam01,
                    lam = lambda[j])
      sig2est[j] <- resj$sig2
      ngest[j] <- length(unique(resj$cluster[,1]))
    }
    BICm <- data.frame(sig2 = sig2est, ng = ngest) %>% 
      mutate(bic = log(sig2) + Cm*ng*(1+ncol(x))*log(nt)/nt)
    
    resm <- SLCC2(indexy = area,y = y_sc, x = cbind(1, x_sc),
                  weights = weights, wtilde = wtilde,betam0 = betam01,
                  lam = lambda[which.min(BICm$bic)])
    
    refit_SR <- estSAE_SR(obj = resm,indexy = area,y = y,
                           x = cbind(1, x), Xbar = cbind(1, Xbar), 
                           wts = wts, N = N, 
                           model ="creg")
  }
  
  if(model == "ccreg")
  {
    ngest = rep(0, nlam)
    
    for(j in 1:nlam)
    {
      resj <- SLCC3(indexy = area,y = y_sc, x = cbind(1,x_sc),group = group,
                    weights = weights, wtilde = wtilde,betam0 = betam01,lam = lambda[j]) 
      sig2est[j] <- resj$sig2
      ngest[j] <- sum(apply(resj$cluster,2, function(x) length(unique(x))))
    }
    BICm <- data.frame(sig2 = sig2est, ng = ngest) %>% 
      mutate(bic = log(sig2) + Cm*ng*log(nt)/nt)
  
    resm <- SLCC3(indexy = area,y = y_sc, x = cbind(1,x_sc),group = group,
                  weights = weights, wtilde = wtilde,betam0 = betam01,
                  lam = lambda[which.min(BICm$bic)]) 
    
    refit_SR <- estSAE_SR(obj = resm,indexy = area,y = y,
                           x = cbind(1, x), Xbar = cbind(1, Xbar), 
                           wts = wts, N = N, 
                           model ="ccreg", group = group)
  }

  out <- list(resm = resm, BICm = BICm,
              refit = refit_SR, cluster = resm$cluster, 
              estY =  refit_SR$estY)
  return(out)
}