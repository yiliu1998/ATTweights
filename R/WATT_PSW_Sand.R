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
