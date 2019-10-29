## This script has the following functions:
##  Generating random numbers:
##   myrexp, myrnorm (N(0,1) using Box-Muller), mymvrnorm (Rejection Sampling),
##    mymvrnorm.

##  Other:  getAIC, getBIC, olscv (will need data preprocessing)

## Exponential

#  F^-1(u) = log(1-u)/-lambda

myrexp <- function(n,mu)
{
  u <- runif(n) # generate a realization of U
  x.list <- log(1 - u) / - mu # generate a realization of X ~ exp(mu)
  return(x.list)
  
}




## The following function uses the Box-Muller
## method to generate a realization of 
## Z_1,..., Z_n iid N(0,1)
## Argument:
##  n, the sample size
##
##  Returns a vector of n entries that has
##  the realization of Z_1,..., Z_n
myrnorm=function(n)
{
  ## is n odd?
  odd=(n %% 2)!=0 
  if(odd) 
  {
    ## add 1 to n 
    n=n+1
  }
  ##  perform the Box-Muller method
  u1=runif(n/2)
  u2=runif(n/2)
  tmp=sqrt(-2*log(u1))
  x.list=c(tmp*cos(2*pi*u2), tmp*sin(2*pi*u2) )
  if(odd)
  {
    ## if n was initially odd, remove an entry
    x.list=x.list[-n]
  }
  return(x.list)
}

## The following function uses rejection sampling 
## to generate a realization of T_1,..., T_n iid with the
## truncated N(mu, sigma^2) distribution on (a,b)
##
##  Arguments:
##    n, the sample size
##    mu, the mean of the un-truncated distribution
##    sigma, the standard deviation of the un-truncated distribution
##    a, the left endpoint
##    b, the right endpoint
##
##  Returns a vector of n entries that has
##  the realization of T_1,..., T_n
myrtnorm = function(n, mu, sigma, a, b)
{
  t.list = numeric(n)
  for(i in 1:n)
  {
    accept=FALSE
    while(!accept)
    {
      ## generate a realization of Z~N(mu,sigma^2)
      z=mu+sigma*rnorm(1)
      accept=(a <= z) & (z <= b) 
    }
    ## at this point, z is a realization of T_i
    t.list[i]=z
  }
  return(t.list)
}


## Multivariate Normal

mymvrnorm <- function(n,mu,Sigma)
{
  p <- length(mu)  # p variate
  # generate z1,..,zn ~ N(0,1) and use these to generate N(mu,Sigma)
  z <- matrix(rnorm(n*p),nrow = n,ncol = p) 
  # get eigen vectors and values
  eig <- eigen(Sigma,symmetric = TRUE)
  # compute Sigma^.5
  sig.sqrt <- eig$vectors %*% diag(eig$values^.5) %*% t(eig$vectors) 
  x <- rep(1,n) %*% t(mu) + z %*%  sig.sqrt
  return(x)
}


## Generate data from G different multivariates for cluster analysis



### AIC/BIC

getAIC = function(x,y){
  
  n = dim(x)[1]
  p = dim(x)[2]
  
  beta.hat=lm.fit(x=X,y=Y)$coefficients
  rss = sum((y - x%*%beta.hat)^2)
  aic = n * log(2*pi) + n * log(rss/n) + n + 2 * (p+1)
  return(aic)
  
}

getBIC = function(x,y){
  
  n = dim(x)[1]
  p = dim(x)[2]
  
  beta.hat=lm.fit(x=X,y=Y)$coefficients
  rss = sum((y - x%*%beta.hat)^2)
  aic = n * log(2*pi) + n * log(rss/n) + n + (p+1)*log(n)
  return(aic)
  
}

olscv = function(x,y,k)
{
  
  # partition the data k times
  folds <- cut(seq(1,nrow(x)),breaks=k,labels=FALSE)
  
  results = numeric(k)
  for(i in 1:k){
    
    testIndexes <- which(folds==i,arr.ind=TRUE)
    
    Xtrain = as.matrix(x[-testIndexes,]) 
    Ytrain = as.matrix(y[-testIndexes])
    
    # fit the linear model
    beta_hat = lm.fit(x=Xtrain,y=Ytrain)$coefficients
    
    # make predictions
    Xtest = as.matrix(x[testIndexes,])
    Ytest = as.matrix(y[testIndexes])
    
    preds = Xtest %*% beta_hat
    
    avg_sq_err = sum((Ytest - preds)^2)/length(Ytest)
    
    results[i] = avg_sq_err
  }
  return(results)
}

## Ridge Binomia Logistic Regression
ilogit=function(X) return(exp(X)/(1 + exp(X)))
ridge.blr=function(X,y,n.list,lam,m=NULL,tol=1e-7,maxit,
                   quiet=FALSE,b.start=NULL){
  
  p = dim(X)[2]
  
  if (is.null(m)) m = c(0,rep(1,p-1))
  
  lam.m = lam*m
  lam.M = diag(lam.m)
  X.t.nlist.y = crossprod(X,n.list*y)
  
  if (is.null(b.start)) {
    b = rep(0,p)
  } else {
    b=b.start
  }
  k=0
  newtdir = tol + 1
  
  while ((k <=maxit) & sum(abs(newtdir)) > tol){
    
    k = k + 1
    pi.t = ilogit(as.numeric(X %*% b))
    W = diag(n.list * pi.t * (1 - pi.t))
    mingrad = X.t.nlist.y - crossprod(X,n.list * pi.t) - lam.m * b
    Hess = crossprod(X,W%*%X) + lam.M
    newtdir = qr.solve(Hess,mingrad)
    b = b + newtdir
    
  }  
  b = as.numeric(b)  
  return(list('b'=b,'total.iterations'=k,'final.gradient' = -mingrad))  
}



### Visualization (ggplot)
library(tidyverse)
bplot_fill = function(data,x,y,fill,xlab = x,ylab=y,...){
  ggplot(data,aes(!!sym(x),!!sym(y),fill= !!sym(fill))) +
    geom_bar(stat='identity',position='dodge') +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    ylab(ylab) + xlab(xlab) +
    #guides(fill=guide_legend(title = 'Treatment')) +
    theme(panel.grid = element_blank()) +
    theme(axis.text=element_text(size=10),axis.title=element_text(size=10), 
          panel.background = element_rect(colour = "black", fill = "white"), 
          strip.text = element_text(size=10)) 
  
}