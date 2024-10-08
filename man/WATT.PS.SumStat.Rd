% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SumStat.R
\name{WATT.PS.SumStat}
\alias{WATT.PS.SumStat}
\title{Summary statistics in propensity score analysis for WATT estimation}
\usage{
WATT.PS.SumStat(
  y,
  z,
  X,
  alpha = 0.05,
  epsilon = 0.001,
  weight = "att",
  trt.SL.library = c("SL.glm"),
  SumStat = "PS.plot"
)
}
\arguments{
\item{y}{the outcome variable, a numeric vector}

\item{z}{the binary treatment assignment, a vector, values need to be 0 or 1}

\item{X}{predictors for the propensity score (confounders), a matrix or a vector (if only one covariate)}

\item{alpha}{trimming or truncation threshold, if using ATT trimming/truncation}

\item{epsilon}{the variance parameter for the smooth trimming method}

\item{weight}{propensity score weighting methods in weighted ATT estimand, valued from
"att", "trimming", "re-est trimming", "smooth trimming", "truncation",
"overlap", "matching", and "entropy";
The "att" means the conventional ATT estimand;
The "trimming" and "re-est trimming": "re-est trimming" will re-estimate the propensity score on the trimmed sample,
where the "trimming" will not;
For "trimming", "re-est trimming" and "truncation", the `alpha` argument must be assigned;
For "smooth trimming", the `epsilon` argument must also be assigned}

\item{trt.SL.library}{methods for fitting the generalized propensity score (chosen from `SuperLearner` package)}

\item{SumStat}{summary statistics need to return, chosen from
"PS.plot" - the propensity score plots by group;
"cov.bal" - covariate balancing by absolute standard mean difference (by the chosen `weight`)}
}
\description{
Summary statistics in propensity score analysis for WATT estimation
}
