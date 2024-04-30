WATT <- function(y, z, X, alpha=0, epsilon=.001, method="att"){
  
  # y: outcome variable
  # z: treatment status (needs to be binary with 0, 1 values)
  # X: predictors for the propensity score
  # alpha: trimming or truncation threshold
  # epsilon: for smooth trimming only, parameter for variance of normal CDF in 
  #          smooth trimming's tilting function
  # method: att, trimming, re-est trimming, smooth trimming, truncation, overlap
  
  # estimate the propensity score via glm
  fit <- glm(z ~ X, family = binomial(link = "logit"))
  e.h <- as.numeric(fit$fitted.values)
  
  y1.h <- sum(y*z)/sum(z)
  if(method!="re-est trimming") {
    if(method=="att")             { tilting <- 1 }
    
    if(method=="overlap")         { tilting <- e.h*(1-e.h) }
    if(method=="matching")        { tilting <- pmin(e.h, 1-e.h) }
    if(method=="entropy")         { tilting <- -e.h*log(e.h) - (1-e.h)*log(1-e.h) }
    
    if(method=="trimming")        { tilting <- I(e.h<=1-alpha) }
    if(method=="smooth trimming") { tilting <- pnorm(1-e.h-alpha, 0, epsilon) }
    if(method=="truncation")      { tilting <- I(e.h<=1-alpha) + (1-alpha)/alpha*(e.h>1-alpha)*(1-e.h)/e.h }
    
    tau <- y1.h - sum(y*(1-z)*e.h/(1-e.h)*tilting)/sum((1-z)*e.h/(1-e.h)*tilting) 
  }
  
  if(method=="re-est trimming") {
    trim <- as.numeric(z==0 & e.h<=1-alpha) + as.numeric(z==1)
    trim.ind <- trim==1
    X <- X[trim.ind,]
    y <- y[trim.ind]
    z <- z[trim.ind]
    
    # re-fit the propensity score model
    fit <- glm(z ~ X, family = binomial(link = "logit"))
    e.h <- as.numeric(fit$fitted.values)
    tau <- sum(y*z)/sum(z) - sum(y*(1-z)*e.h/(1-e.h))/sum((1-z)*e.h/(1-e.h)) 
  }
  return(tau)
}
