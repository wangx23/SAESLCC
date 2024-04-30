####### get area estimates based on regression models ####
# obj is the object of SLCC
# Xbar population mean of X for each area
# Zbar population mean of Z for each area
# wts sampling weights
# N population area size
### based on the cluster information in obj to refit the model 
### and then construct the estimator 

### model is the model used to construct estimator
## three options: "intercept" (model 1), "creg" (model 2), "ccreg" (model 3)
# group is coordinate group for model 3

estSAE_SR <- function(obj, indexy, y, x, Xbar, z=NULL, Zbar = NULL, 
                      wts, N, model, group=NULL)
{
  cluster_est <- obj$cluster
  nivec <- as.numeric(table(indexy))
  Nbar <- mean(N)
  wtilde  <- rep(1/Nbar, nivec) * wts ### adjusted weight in the algorithm
  if(model =="intercept")
  {
    res_refit <- refit_m1(indexy = indexy, y = y, x=x, z=z,
                          cluster = cluster_est, wtilde = wtilde)
  }
  if(model=="creg")
  {
    res_refit <- refit_m2(indexy = indexy, y=y, x=x, cluster=cluster_est, wtilde)
  }
  if(model =="ccreg")
  {
    res_refit <- refit_m3(indexy=indexy, y=y, x=x, group=group, 
                          clustermat = cluster_est,wtilde)
  }
  
  uindexy <- unique(indexy)
  yhat <- res_refit$yhat
  
  subfun <- function(ii)
  {
    indexii <- indexy == uindexy[ii]
    yii <- y[indexii]
    yhatii <- yhat[indexii]
    
    return(sum((yii - yhatii)*wts[indexii]))
  }
  estp1 <- unlist(lapply(1:length(uindexy),subfun))
  
  eta_est <- res_refit$eta
  beta_est <- res_refit$beta
  
  if(is.null(Zbar))
  {
    estY <- rowSums(Xbar * beta_est) + estp1/N
  }else{
    estY <- Zbar %*% eta_est + rowSums(Xbar * beta_est) + estp1/N
  }
  
  out <- list(estY = estY, refit = res_refit, model = model)
  return(out)
}