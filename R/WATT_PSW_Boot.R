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
  X <- as.data.frame(X)
  fit <- SuperLearner(Y=z, X=X, SL.library=trt.SL.library, family=binomial())
  e.h <- predict(fit, X, type="response")$pred

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

