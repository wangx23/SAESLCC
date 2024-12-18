# SAESLCC
library(devtools)

install_github("wangx23/SAESLCC")


#### run a simulated example 



```
data(simdata)

library(SAESLCC)
area <- dat_sample$area
y_sample <- dat_sample$y
x_sample <- as.matrix(dat_sample$x1)
wts <- dat_sample$wts

plot(x_sample, y_sample)

##### needs to specify lambda values 

lamvec1 <- seq(0.02,0.4,length.out = 50)
lamvec2 <- seq(0.02,0.5,length.out = 50)
lamvec3 <- seq(0.02,0.3,length.out = 50)


res1 <- SLCC_SR(area = area, y = y_sample, x = x_sample,Xbar = Xbar,
                wts = wts, N = N, model = "intercept",lambda = lamvec1)
                


res2 <- SLCC_SR(area = area, y = y_sample, x = x_sample,Xbar = Xbar,
                wts = wts, N = N, model = "creg",lambda = lamvec2)
                

plot(res2$refit$refit$beta) ### estimates
res3 <- SLCC_SR(area = area, y = y_sample, x = x_sample,Xbar = Xbar,
                wts = wts, N = N, model = "ccreg",lambda = lamvec3, 
                group = c(1,2),standardize =TRUE, center = FALSE)


cbind(res1$estY, res2$estY, res3$estY) #estimates 


```


