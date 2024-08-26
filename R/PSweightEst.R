#' Point estimation of a chosen WATT estimands, using propensity score weighting estimator
#' @param y the outcome variable, a numeric vector
#' @param z the binary treatment assignment, a vector, values need to be 0 or 1
#' @param X predictors (confounders) for the propensity score, a matrix or a vector (if only one covariate)
#' @param alpha trimming or truncation threshold, if using ATT trimming/truncation
#' @param epsilon the variance parameter for the smooth trimming method
#' @param weight propensity score weighting methods in weighted ATT estimand, valued from
#'               "att", "trimming", "re-est trimming", "smooth trimming", "truncation",
#'               "overlap", "matching", and "entropy";
#'               The "att" means the conventional ATT estimand;
#'               The "trimming" and "re-est trimming": "re-est trimming" will re-estimate the propensity score on the trimmed sample,
#'               where the "trimming" will not;
#'               For "trimming", "re-est trimming" and "truncation", the `alpha` argument must be assigned;
#'               For "smooth trimming", the `epsilon` argument must also be assigned
#' @param trt.SL.library methods for fitting the propensity score (chosen from `SuperLearner` package)
WATT.PSW <- function(y, z, X,
                     alpha=0.05, epsilon=0.001, weight="att",
                     trt.SL.library="SL.glm"){

  library(SuperLearner)
  library(dplyr)
  library(ggplot2)

  # estimate the propensity score
  if(trt.SL.library=="SL.glm") {
    fit <- glm(z ~ X, family = binomial(link = "logit"))
    e.h <- as.numeric(fit$fitted.values, type="response")
  } else {
    X <- as.data.frame(X)
    fit <- SuperLearner(Y=z, X=X, SL.library=trt.SL.library, family=binomial())
    e.h <- predict(fit, X, type="response")$pred
  }

  y1.h <- sum(y*z)/sum(z)

  if(weight!="re-est trimming") {
    if(weight=="att")             { tilting <- 1 }

    if(weight=="overlap")         { tilting <- e.h*(1-e.h) }
    if(weight=="matching")        { tilting <- pmin(e.h, 1-e.h) }
    if(weight=="entropy")         { tilting <- -e.h*log(e.h) - (1-e.h)*log(1-e.h) }

    if(weight=="trimming")        { tilting <- I(e.h<=1-alpha) }
    if(weight=="smooth trimming") { tilting <- pnorm(1-e.h-alpha, 0, epsilon) }
    if(weight=="truncation")      { tilting <- I(e.h<=1-alpha) + (1-alpha)/alpha*I(e.h>1-alpha)*(1-e.h)/e.h }

    tau <- y1.h - sum(y*(1-z)*e.h/(1-e.h)*tilting)/sum((1-z)*e.h/(1-e.h)*tilting)
    w1 <- z
    w0 <- (1-z)*tilting*e.h/(1-e.h)
  }

  if(weight=="re-est trimming") {
    trim <- as.numeric(z==0 & e.h<=1-alpha) + as.numeric(z==1)
    trim.ind <- trim==1
    if(class(X)[1]!="matrix") {
      X <- X[trim.ind]
    } else {
      X <- X[trim.ind, ] }
    y <- y[trim.ind]
    z <- z[trim.ind]

    # re-fit the propensity score
    fit <- SuperLearner(Y=z, X=X, SL.library=trt.SL.library, family=binomial())
    e.h <- predict(fit, X)$pred
    tau <- sum(y*z)/sum(z) - sum(y*(1-z)*e.h/(1-e.h))/sum((1-z)*e.h/(1-e.h))
    w1 <- z
    w0 <- (1-z)*e.h/(1-e.h)
  }
  return(list(tau=tau, w1=w1, w0=w0))
}

#' Bootstrap variance estimation for WATT weighting estimator
#' @param y the outcome variable, a numeric vector
#' @param z the binary treatment assignment, a vector, values need to be 0 or 1
#' @param X predictors for the propensity score (confounders), a matrix or a vector (if only one covariate)
#' @param alpha trimming or truncation threshold, if using ATT trimming/truncation
#' @param epsilon the variance parameter for the smooth trimming method
#' @param weight propensity score weighting methods in weighted ATT estimand, valued from
#'               "att", "trimming", "re-est trimming", "smooth trimming", "truncation",
#'               "overlap", "matching", and "entropy";
#'               The "att" means the conventional ATT estimand;
#'               The "trimming" and "re-est trimming": "re-est trimming" will re-estimate the propensity score on the trimmed sample,
#'               where the "trimming" will not;
#'               For "trimming", "re-est trimming" and "truncation", the `alpha` argument must be assigned;
#'               For "smooth trimming", the `epsilon` argument must also be assigned
#' @param trt.SL.library methods for fitting the generalized propensity score (chosen from `SuperLearner` package)
#' @param N.boot number of bootstrap replicates
#' @param seed seed for random numbers, used for bootstrap sampling
#' @param conf.level level of constructed confidence interval (using normal approximation); default is 0.95
WATT.PSW.bootstrap <- function(y, z, X,
                               alpha=0.05, epsilon=0.001, weight="att",
                               trt.SL.library=c("SL.glm"),
                               N.boot=200,
                               seed=0509,
                               conf.level=0.95){

  tau <- WATT.PSW(y=y, z=z, X=X,
                  alpha=alpha, epsilon=epsilon, weight=weight, trt.SL.library=trt.SL.library)$tau
  set.seed(seed=seed)
  seeds <- round(runif(N.boot, min=0, max=1e5*N.boot))
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
    tau.boot[i] <- WATT.PSW(y=y[boot], z=z[boot], X=X.boot,
                            alpha=alpha, epsilon=epsilon, weight=weight, trt.SL.library=trt.SL.library)$tau

    if(i%%50==0) print(paste0("bootstrap ", i, " is done."))
  }
  sd.tau <- sd(tau.boot)
  quant <- qnorm(1-(1-conf.level)/2)
  lwr.tau <- tau - sd.tau*quant
  upr.tau <- tau + sd.tau*quant
  df <- data.frame(tau=tau, sd=sd.tau, lwr=lwr.tau, upr=upr.tau)
  return(df)
}

#' Sandwich variance estimation for WATT weighting estimator (only allows logistic regression model for PS)
#' @param y the outcome variable, a numeric vector
#' @param z the binary treatment assignment, a vector, values need to be 0 or 1
#' @param X predictors for the propensity score (confounders), a matrix or a vector (if only one covariate)
#' @param alpha trimming or truncation threshold, if using ATT trimming/truncation
#' @param epsilon the variance parameter for the smooth trimming method
#' @param weight propensity score weighting methods in weighted ATT estimand, valued from
#'               "att", "trimming", "smooth trimming", "truncation", "overlap", "matching", and "entropy";
#'               The "att" means the conventional ATT estimand;
#'               For "trimming" and "truncation", the `alpha` argument must be assigned; default is 0.05;
#'               For "smooth trimming", the `epsilon` argument must also be assigned; default is 0.001
#' @param conf.level level of constructed confidence interval (using normal approximation); default is 0.95
WATT.PSW.sandwich <- function(y, z, X,
                              alpha=0.05, epsilon=0.001, weight="att", conf.level=0.95){

  fit <- glm(z ~ X, family = binomial(link = "logit"))
  e.h <- as.numeric(fit$fitted.values)
  n <- length(z)
  V <- cbind(1, X)

  if(weight=="att") {
    g.h <- 1
    dg.h <- 0
  }
  if(weight=="overlap") {
    g.h <- e.h*(1-e.h)
    dg.h <- (1-2*e.h)*e.h*(1-e.h) * V
  }
  if(weight=="matching") {
    g.h <- pmin(e.h, 1-e.h)
    dg.h <- (I(e.h<0.5)-I(e.h>0.5))*e.h*(1-e.h) * V
  }
  if(weight=="entropy") {
    g.h <- -e.h*log(e.h) - (1-e.h)*log(1-e.h)
    dg.h <- (log(1-e.h)-log(e.h))*e.h*(1-e.h) * V
  }
  if(weight=="trimming") {
    g.h <- I(e.h<=1-alpha)
    dg.h <- 0
  }
  if(weight=="smooth trimming") {
    g.h <- pnorm(1-e.h-alpha, 0, epsilon)
  }
  if(weight=="truncation") {
    g.h <- I(e.h<=1-alpha) + (1-alpha)/alpha*I(e.h>1-alpha)*(1-e.h)/e.h
    dg.h <- (1-alpha)/alpha*I(e.h>1-alpha)*(e.h-1)/e.h * V
  }

  omega0.h <- g.h*e.h/(1-e.h)
  mu1.h <- z*y / sum(z)
  mu0.h <- (1-z)*y*omega0.h / sum((1-z)*omega0.h)
  tau <- sum(mu1.h)-sum(mu0.h)

  A.11 <- crossprod(sqrt(e.h*(1-e.h)) * V) / n
  A.22 <- mean(z)
  A.31 <- apply((1-z)*(-e.h/(1-e.h))*(dg.h + g.h*V)*(y-mu0.h), 2, mean)
  A.33 <- mean( omega0.h*(1-z) )

  C.11 <- A.11
  C.21 <- rbind(0, A.31)
  C.22 <- diag(c(A.22, A.33))
  cA <- c(1,-1) %*% cbind(-solve(C.22)%*%C.21%*%solve(C.11), solve(C.22))

  psi <- cbind((z-e.h)*V,  z*(y-mu1.h),  (1-z)*omega0.h*y-e.h*g.h*mu0.h )
  B <- crossprod(psi) / n
  VAR <- cA %*% B %*% t(cA) / n
  sd.tau <- as.numeric(sqrt(VAR))

  quant <- qnorm(1-(1-conf.level)/2)
  lwr.tau <- tau - sd.tau*quant/sqrt(n)
  upr.tau <- tau + sd.tau*quant/sqrt(n)
  df <- data.frame(tau=tau, sd=sd.tau, lwr=lwr.tau, upr=upr.tau)
  return(df)
}
