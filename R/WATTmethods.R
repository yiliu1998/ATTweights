WATT <- function(y, z, X, 
                 alpha=0, epsilon=.001, 
                 method="att"){
  
  # y: outcome variable, a numeric vector
  # z: treatment status (needs to be binary with 0, 1 values), a vector
  # X: predictors for the propensity score, a matrix or a vector (if only one covariate)
  # alpha: trimming or truncation threshold
  # epsilon: for smooth trimming only, parameter for variance of normal CDF in 
  #          smooth trimming's tilting function
  # method: att, trimming, re-est trimming, smooth trimming, truncation, overlap, matching, entropy
  
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
    if(method=="truncation")      { tilting <- I(e.h<=1-alpha) + (1-alpha)/alpha*I(e.h>1-alpha)*(1-e.h)/e.h }
    
    tau <- y1.h - sum(y*(1-z)*e.h/(1-e.h)*tilting)/sum((1-z)*e.h/(1-e.h)*tilting) 
    w1 <- z
    w0 <- (1-z)*tilting*e.h/(1-e.h)
  }
  
  if(method=="re-est trimming") {
    trim <- as.numeric(z==0 & e.h<=1-alpha) + as.numeric(z==1)
    trim.ind <- trim==1
    if(class(X)[1]!="matrix") { 
      X <- X[trim.ind] 
    } else { 
      X <- X[trim.ind, ] }
    y <- y[trim.ind]
    z <- z[trim.ind]
    
    # re-fit the propensity score model
    fit <- glm(z ~ X, family = binomial(link = "logit"))
    e.h <- as.numeric(fit$fitted.values)
    tau <- sum(y*z)/sum(z) - sum(y*(1-z)*e.h/(1-e.h))/sum((1-z)*e.h/(1-e.h)) 
    w1 <- z
    w0 <- (1-z)*e.h/(1-e.h)
  }
  return(list(tau=tau, w1=w1, w0=w0))
}

WATT.bootstrap <- function(y, z, X, 
                           alpha=0, epsilon=.001, 
                           method="att",
                           N.boot=200, 
                           seed=0509,
                           conf.level=0.05){
  tau <- WATT(y, z, X, alpha, epsilon, method)$tau
  set.seed(seed=seed)
  seeds <- round(runif(N.boot*1e2, min=0, max=1e5*N.boot))
  n <- length(z)
  tau.boot <- c()
  for(i in 1:N.boot) {
    set.seed(seeds[i])
    boot <- sample(1:n, n, replace=TRUE)
    if(class(X)[1]!="matrix") {
      X.boot <- X[boot]
    } else {
      X.boot <- X[boot,]
    }
    tau.boot[i] <- WATT(y=y[boot], z=z[boot], X=X.boot, 
                        alpha, epsilon, method)$tau
  }
  sd.tau <- sd(tau.boot)
  quant <- qnorm(1-conf.level/2)
  lwr.tau <- tau - sd.tau*quant/sqrt(n)
  upr.tau <- tau + sd.tau*quant/sqrt(n)
  df <- data.frame(tau=tau, sd=sd.tau, lwr=lwr.tau, upr=upr.tau)
  return(df)
}

WATT.SumStat <- function(y, z, X, 
                         alpha=0, epsilon=.001, 
                         method="att",
                         SumStat="PS.plot") {
  
  # SumStat: desired summary statistics, choosing from
  # --- PS.plot - the propensity score plots by group
  # --- cov.bal - covariate balancing by absolute standard mean difference 
  #               using the chosen weighting method
  # --- all - return all summary statistics above
  
  library(ggplot2)
  
  # estimate the propensity score via glm
  fit <- glm(z ~ X, family = binomial(link = "logit"))
  e.h <- as.numeric(fit$fitted.values)
  df <- data.frame(e.h=e.h, z=z)
  ps.plot <- ggplot(df, aes(x=e.h, fill=factor(z))) + 
    geom_histogram(position="identity", alpha=0.5, color="black", bins=35) + 
    labs(x="Estimated propensity score", y="Count") +
    scale_fill_manual(values=c("darkblue", "gray"), name="Group", labels=c("Control","Treated")) + 
    theme_minimal()
  
  # balancing checking by the chosen weighting method
  w1 <- WATT(y, z, X, alpha, epsilon, method)$w1
  w0 <- WATT(y, z, X, alpha, epsilon, method)$w0
  
  if(class(X)[1]!="matrix") {
    x1.mean <- sum(w1*X)/sum(w1)
    x0.mean <- sum(w0*X)/sum(w0)
    s1.sq <- sum(w1)/((sum(w1))^2 - sum(w1^2)) * sum(w1*(X-x1.mean)^2)
    s0.sq <- sum(w0)/((sum(w0))^2 - sum(w0^2)) * sum(w0*(X-x0.mean)^2)
    ASMD <- abs((x1.mean-x0.mean)/sqrt(s1.sq/2+s0.sq/2))
  } else {
    ASMD <- c()
    for(j in 1:ncol(X)) {
      x1.mean <- sum(w1*X[,j])/sum(w1)
      x0.mean <- sum(w0*X[,j])/sum(w0)
      s1.sq <- sum(w1)/((sum(w1))^2 - sum(w1^2)) * sum(w1*(X[,j]-x1.mean)^2)
      s0.sq <- sum(w0)/((sum(w0))^2 - sum(w0^2)) * sum(w0*(X[,j]-x0.mean)^2)
      ASMD[j] <- abs((x1.mean-x0.mean)/sqrt(s1.sq/2+s0.sq/2))
    }
    if(is.null(colnames(X))) { 
      X.names <- paste0("X", 1:ncol(X)) 
      } else {
      X.names <- colnames(X)
      }
    ASMD <- data.frame(matrix(ASMD, nrow=1))
    colnames(ASMD) <- X.names
  }
  if(SumStat=="PS.plot") { return(ps.plot) }
  if(SumStat=="cov.bal") { return(ASMD) }
  if(SumStat=="all") { return(list(ps.plot=ps.plot, ASMD=ASMD)) }
}
