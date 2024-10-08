# ATTweights
`ATTweights` is an R package with illustrative examples for weighted ATT methods (a generalization to the conventional ATT estimand). The main methodology paper of `ATTweights` is [Liu et al. (2024)](https://journals.sagepub.com/doi/10.1177/09622802241269646). The current verison of this package can implment the following estimators: 

* Propensity score (PS) weighting estimator for WATT, with (i) bootstrap; (ii) sandwich variance estimations. For (i), we allow users to specify different PS and outcome regression (OR) models used in `SuperLearner` R package (see <a href="https://academic.oup.com/ectj/article/21/1/C1/5056401">Chernozhukov et al. (2018)</a>). For (ii), we only allow linear model for the OR models and logistic regression for the PS model. Usually, the sandwich variance estimation is designed for parametric models. When using machine learning for PS and OR models, one can also consider a direct variance estimation by the mean square of influence functions, as the bootstrap can be time-consuming. However, for weighting estimator, its influence function under different models is still not clear in literature, so we have not yet included this variance option in the package. 

* Augmented estimator for WATT with only sandwich variance estimator, allowing only linear models for OR and logistic regression for PS.

<strong>Remark for sandwich variance:</strong> The sandwich variance formula relies on the parametric models used for PS and OR, and in our current package version, we have not extended it to other models rather than logistic regression for PS and linear regreesion for OR models (e.g., when the outcome is binary and one wishes to use logistic regression to model it, our sandwich variance is not applicable). We will think about extending it to other models, but currently we suggest the use of bootstrap variance estimation when the sandwich variance is not applicable in your case. 

## Installation
To install the latest version of the R package from GitHub, please run following commands in R:

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("yiliu1998/ATTweights")
```

## Use

Let us generate the following data to test our method. You can copy and paste the code on your local R studio. 

```r
n <- 2000
X1 <- rnorm(n, mean=1, sd=2)
X2 <- rnorm(n, mean=0, sd=3)
X3 <- rt(n, df=10)
X4 <- rbinom(n, size=20, prob=0.6*abs(sin(X1)))
X5 <- X1*X2
X6 <- X3^2
X7 <- X4*X3

expit <- function(x) 1/(1+exp(-x))
beta <- c(-0.5, rep(0.1, 7))
X <- cbind(X1,X2,X3,X4,X5,X6,X7)
X.ps <- X.out <- cbind(1,X)

ps <- expit(X.ps%*%beta)
Z <- rbinom(n, 1, ps)
mean(Z)

alpha <- c(5, rep(1,7))
Y <- X.out%*%alpha + Z*15 + rnorm(n)
```

Running the following code, we can get a propensity score distribution plot by treatment group. 

```r
library(ATTweights) # assuming the latest package is installed correctly
WATT.PS.SumStat(y=Y, z=Z, X=X)
```

Running the following two code, we can get the point estimates, standard errors and confidence intervals by bootstrap (500 replicates) of OWATT and ATT, respectively. This code uses PS weighting estimator.  

```r
WATT.PSW.bootstrap(y=Y, z=Z, X=X, weight="overlap", N.boot=500)
WATT.PSW.bootstrap(y=Y, z=Z, X=X, weight="att", N.boot=500)
```

## Contact 
The R code is maintained by Yi Liu (Please feel free to reach out at  yi.liu.biostat@gmail.com, if you have any questions). 

## Reference
Please cite the following paper:

Liu Y, Li H, Zhou Y, Matsouaka RA. Average treatment effect on the treated, under lack of positivity. Statistical Methods in Medical Research. 2024;0(0). [doi:10.1177/09622802241269646](https://journals.sagepub.com/doi/10.1177/09622802241269646)
