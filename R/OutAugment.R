#' Augmented WATT estimation and sandwich variance (only linear model for outcomes and logistic regression for PS allowed)
#' @param y the outcome variable, a numeric vector
#' @param z the binary treatment assignment, a vector, values need to be 0 or 1
#' @param X predictors (confounders) for the propensity score, a matrix or a vector (if only one covariate)
#' @param X.out predictors (confounders) for the two potential outcomes, a matrix or a vector (if only one covariate)
#' @param alpha trimming or truncation threshold, if using ATT trimming/truncation
#' @param epsilon the variance parameter for the smooth trimming method
#' @param weight propensity score weighting methods in weighted ATT estimand, valued from
#'               "att", "trimming", "smooth trimming", "truncation", "overlap", "matching", and "entropy";
#'               The "att" means the conventional ATT estimand;
#'               For "trimming" and "truncation", the `alpha` argument must be assigned; default is 0.05;
#'               For "smooth trimming", the `epsilon` argument must also be assigned; default is 0.001
#' @param conf.level level of constructed confidence interval (using normal approximation); default is 0.95
WATT.Aug.sandwich <- function(y, z, X, X.out,
                     alpha=0.05, epsilon=0.001, weight="att"){

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
    dg.h <- dnorm(1-e.h-alpha, 0, epsilon)
  }
  if(weight=="truncation") {
    g.h <- I(e.h<=1-alpha) + (1-alpha)/alpha*I(e.h>1-alpha)*(1-e.h)/e.h
    dg.h <- (1-alpha)/alpha*I(e.h>1-alpha)*(e.h-1)/e.h * V
  }

  out.ctrl <- lm(y~., data=data.frame(y, X.out)[z==0,])
  out.trt  <- lm(y~., data=data.frame(y, X.out)[z==1,])
  m0.h <- predict(out.ctrl, as.data.frame(X.out))
  m1.h <- predict(out.trt, as.data.frame(X.out))
  omega0.h <- g.h*e.h/(1-e.h)

  mu1.h <- sum( z*(y-m1.h) ) / sum(z)
  mu0.h <- sum( (1-z)*omega0.h*(y-m0.h) ) / sum( (1-z)*omega0.h )
  eta1.h <- sum( e.h*m1.h ) / sum(e.h)
  eta0.h <- sum( e.h*g.h*m0.h ) / sum(e.h*g.h)

  W <- cbind(1, X.out)

  A.11 <- crossprod(sqrt(e.h*(1-e.h)) * V) / n
  A.22 <- crossprod(z * W) / n
  A.33 <- crossprod((1-z) * W) / n
  A.41 <- apply(-e.h*(1-e.h)*(m1.h-eta1.h) * V, 2, mean)
  A.42 <- apply(-e.h * W, 2, mean)
  A.44 <- mean(e.h)
  A.51 <- apply(-e.h*((1-e.h)*g.h * V + dg.h)*(m0.h-eta0.h), 2, mean)
  A.53 <- apply(-e.h*g.h * W, 2, mean)
  A.55 <- mean(e.h*g.h)
  A.62 <- apply(z * W, 2, mean)
  A.66 <- mean(z)
  A.71 <- apply(-(1-z)*e.h/(1-e.h)*(dg.h + g.h * V)*(y-m0.h-mu0.h), 2, mean)
  A.73 <- apply((1-z)*omega0.h * W, 2, mean)
  A.77 <- mean( (1-z)*omega0.h )

  C.11 <- rbind(cbind(A.11, matrix(0,nrow(A.11),(ncol(A.22)+ncol(A.33)))),
                cbind(matrix(0,nrow(A.22),ncol(A.11)), A.22, matrix(0,nrow(A.22),ncol(A.33))),
                cbind(matrix(0,nrow(A.33),(ncol(A.11)+ncol(A.22))), A.33)   )
  C.21 <- rbind(c(A.41, A.42, rep(0, ncol(A.33))),
                c(A.51, rep(0, ncol(A.22)), A.53),
                c(rep(0, ncol(A.11)), A.62, rep(0, ncol(A.33))),
                c(A.71, rep(0, ncol(A.22)), A.73)   )
  C.22 <- diag(c(A.44, A.55, A.66, A.77))

  cA <- c(1,-1,1,-1) %*% cbind(-solve(C.22)%*%C.21%*%solve(C.11), solve(C.22))

  psi <- cbind((z-e.h)*V,
               c( z * (y-m1.h) )* W,
               c( (1-z) * (y-m0.h) )* W,
               e.h*(m1.h-eta1.h),
               e.h*g.h*(m0.h-eta0.h),
               z*(y-m1.h-mu1.h),
               (1-z)*omega0.h*(y-m0.h)-e.h*g.h*mu0.h )
  B <- crossprod(psi) / n

  VAR <- cA %*% B %*% t(cA) / n
  sd.tau <- as.numeric(sqrt(VAR))

  quant <- qnorm(1-(1-conf.level)/2)
  lwr.tau <- tau - sd.tau*quant/sqrt(n)
  upr.tau <- tau + sd.tau*quant/sqrt(n)
  df <- data.frame(tau=tau, sd=sd.tau, lwr=lwr.tau, upr=upr.tau)
  return(df)
}
