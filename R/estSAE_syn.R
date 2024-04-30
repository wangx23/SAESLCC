####### get area estimates based on synthetic regression models ####
# obj is the object of SLCC
# Xbar population mean of X for each area
# Zbar population mean of Z for each area
# wts sampling weights
# N population area size
estSAE_syn <- function(obj, indexy, y, x, Xbar, z=NULL, Zbar = NULL, 
                       wts, N, model, group=NULL)
{
  cluster_est <- obj$cluster
  nivec <- as.numeric(table(indexy))
  Nbar <- rep(mean(N), length(N))
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

  eta_est <- res_refit$eta
  beta_est <- res_refit$beta
  
  if(is.null(Zbar))
  {
    estY <- rowSums(Xbar * beta_est)
  }else{
    estY <- Zbar %*% eta_est + rowSums(Xbar * beta_est)
  }
  out <- list(estY = estY, refit = res_refit, model = model)
  return(out)
}