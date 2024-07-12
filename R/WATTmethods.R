#' Load necessary packages
library(SuperLearner)
library(dplyr)
library(ggplot2)
# listWrappers()

#' Function WATT(): implements the point estimation of a weighted ATT method
#' @param y the outcome variable, a numeric vector
#' @param z the binary treatment assignment, a vector, values need to be 0 or 1
#' @param X predictors for the propensity score (confounders), a matrix or a vector (if only one covariate)
#' @param alpha trimming or truncation threshold, if using ATT trimming/truncation
#' @param epsilon the variance parameter for the smooth trimming method
#' @param weight propensity score weighting methods in weighted ATT estimand, valued from
#'               "att", "trimming", "re-est trimming", "smooth trimming", "truncation",
#'               "overlap", "matching", and "entropy"
#'               The "att" means the conventional ATT estimand
#'               The "trimming" and "re-est trimming": "re-est trimming" will re-estimate the propensity score on the trimmed sample,
#'               where the "trimming" will not
#'               For "trimming", "re-est trimming" and "truncation", the alpha argument must be assigned
#'               For "smooth trimming", the epsilon argument must also be assigned
#' @param trt.SL.library methods for fitting the generalized propensity score (chosen from SuperLearner package)
WATT <- function(y, z, X,
                 alpha=0, epsilon=.001, weight="att",
                 trt.SL.library=c("SL.glm")){

  # estimate the propensity score
  fit <- SuperLearner(Y=z, X=X, SL.library=trt.SL.library, family=binomial())
  e.h <- predict(fit, X)$pred

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

#' Function WATT.bootstrap(): implements the bootstrap variance estimation of the weighting estimator of a WATT estimand
#' @param N.boot number of bootstrap replicates
#' @param seed seed for random numbers, used for bootstrap sampling
#' @param conf.level level of constructed confidence interval (using normal approximation), with default 0.95
#' Other parameters are the same in the WATT() function above
WATT.bootstrap <- function(y, z, X,
                           alpha=0, epsilon=.001, weight="att",
                           trt.SL.library=c("SL.glm"),
                           N.boot=200,
                           seed=0509,
                           conf.level=0.95){

  tau <- WATT(y=y, z=z, X=X,
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
    tau.boot[i] <- WATT(y=y[boot], z=z[boot], X=X.boot,
                        alpha=alpha, epsilon=epsilon, weight=weight, trt.SL.library=trt.SL.library)$tau
  }
  sd.tau <- sd(tau.boot)
  quant <- qnorm(1-(1-conf.level)/2)
  lwr.tau <- tau - sd.tau*quant/sqrt(n)
  upr.tau <- tau + sd.tau*quant/sqrt(n)
  df <- data.frame(tau=tau, sd=sd.tau, lwr=lwr.tau, upr=upr.tau)
  return(df)
}

#' Function WATT.SumStat(): summary statistics in propensity score analysis
#' @param SumStat summary statistics need to return, chosen from
#'                "PS.plot" - the propensity score plots by group
#'                "cov.bal" - covariate balancing by absolute standard mean difference (by the chosen "weight")
#'                "all" - return all summary statistics above
#' Other parameters are the same in the WATT() function above
WATT.SumStat <- function(y, z, X,
                         alpha=0, epsilon=.001,
                         weight="att",
                         trt.SL.library=c("SL.glm"),
                         SumStat="PS.plot") {

 # estimate the propensity score
  fit <- SuperLearner(Y=z, X=X, SL.library=trt.SL.library, family=binomial())
  e.h <- predict(fit, X)$pred

  df <- data.frame(e.h=e.h, z=z)
  ps.plot <- ggplot(df, aes(x=e.h, fill=factor(z))) +
    geom_histogram(position="identity", alpha=0.5, color="black", bins=35) +
    labs(x="Estimated propensity score", y="Count") +
    scale_fill_manual(values=c("darkblue", "gray"), name="Group", labels=c("Control","Treated")) +
    theme_minimal()

  # balancing checking by the chosen weighting method
  w1 <- WATT(y, z, X, alpha, epsilon, weight)$w1
  w0 <- WATT(y, z, X, alpha, epsilon, weight)$w0

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
