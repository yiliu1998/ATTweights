#' Summary statistics in propensity score analysis for WATT estimation
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
#' @param SumStat summary statistics need to return, chosen from
#'                "PS.plot" - the propensity score plots by group;
#'                "cov.bal" - covariate balancing by absolute standard mean difference (by the chosen `weight`)
WATT.PS.SumStat <- function(y, z, X,
                            alpha=0.05, epsilon=0.001,
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
  w1 <- WATT.PSW(y, z, X, alpha, epsilon, weight)$w1
  w0 <- WATT.PSW(y, z, X, alpha, epsilon, weight)$w0

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
