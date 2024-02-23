##### refit based on given cluster information and sampling weights ####

#### three models ###

#### model 1#####
#  z common regression 
# x cluster
# cluster is prespecified cluster structure
#### 
refit_m1<- function(indexy, y, z, x, cluster, wtilde)
{
  ncx <- ncol(x)
  nr <- length(cluster)
  Xm <- matrix(0, nrow(x), nr*ncx)
  uniqxy <- unique(indexy)
  ng <- length(unique(cluster))
  
  for(i in 1:nr)
  {
    Xm[indexy == uniqxy[i],(ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniqxy[i],]
  }
  
  W <- matrix(0, nr, ng)
  W[cbind(1:nr,cluster)] <- 1
  W <- W %x% diag(1,ncx)
  Ux <- cbind(Xm%*%W, z)
  Ux1 <- sqrt(wtilde) *Ux
  y1 <- sqrt(wtilde) *y
  est <- solve(t(Ux1)%*%Ux1)%*%t(Ux1)%*%y1
  
  sig2 <- sum((y1 - Ux1%*%est)^2)/nr
  
  eta_est <- est[-(1:(ng*ncx))]
  beta_est <- matrix(est[1:(ng*ncx)],ng, byrow = TRUE)
  beta_est <- beta_est[cluster,]
  yhat <- Ux %*% est
  
  out <- list(eta = eta_est, beta = beta_est, yhat = yhat, sig2 = sig2)
  return(out)
}


#### model 2 ####
## all x have groups, no common group 
refit_m2 <- function(indexy, y, x, cluster, wtilde)
{
  ncx <- ncol(x)
  nr <- length(cluster)
  Xm <- matrix(0, nrow(x), nr*ncx)
  uniqxy <- unique(indexy)
  ng <- length(unique(cluster))
  
  for(i in 1:nr)
  {
    Xm[indexy == uniqxy[i],(ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniqxy[i],]
  }
  
  W <- matrix(0, nr, ng)
  W[cbind(1:nr,cluster)] <- 1
  W <- W %x% diag(1,ncx)
  Ux <- Xm%*%W
  Ux1 <- sqrt(wtilde) *Ux
  y1 <- sqrt(wtilde) *y
  est <- solve(t(Ux1)%*%Ux1)%*%t(Ux1)%*%y1
  
  sig2 <- sum((y1 - Ux1%*%est)^2)/nr
  
  beta_est <- matrix(est[1:(ng*ncx)],ng, byrow = TRUE)
  beta_est <- beta_est[cluster,]
  yhat <- Ux %*% est
  
  out <- list(beta = beta_est, yhat = yhat, sig2= sig2)
  
  return(out)

}


### model 3 ###
## each coordinate has its own group information
# groupmat is the true group structure
# group is the group index for coordinate 
refit_m3 <- function(indexy, y, x, group, clustermat,wtilde)
{
  ncx <- ncol(x)
  N <- length(y)
  nr <- nrow(clustermat)
  clustermat <- clustermat[,group]
  ngest <- apply(clustermat, 2, function(x){length(unique(x))})
  ngtotal <- sum(ngest)
  
  Xm <- matrix(0, N, nr*ncx)
  uniqxy <- unique(indexy)
  
  for(i in 1:nr)
  {
    Xm[indexy == uniqxy[i],(ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniqxy[i],]
  }
  
  ngestcum <- c(0,cumsum(ngest)[-ncx])
  wfun <- function(i)
  {
    W <- matrix(0, ncx, ngtotal)
    W[cbind(1:ncx,clustermat[i,] + ngestcum)] <- 1
    return(W)
  }
  Wmat <- do.call("rbind",lapply(1:nr, wfun))
  
  Ux <- Xm %*% Wmat
  Ux1 <- sqrt(wtilde) *Ux
  y1 <- sqrt(wtilde) *y
  est <- solve(t(Ux1)%*%Ux1)%*%t(Ux1)%*%y1
  
  sig2 <- sum((y1 - Ux1%*%est)^2)/nr
  yhat <- Ux %*% est
 
  estm <- Wmat %*% est
  beta <- matrix(estm, nr, byrow=TRUE)
  
  
  out <- list(beta=beta, yhat = yhat, sig2 = sig2)

  return(out)
}



