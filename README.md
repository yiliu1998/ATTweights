# ATTweights
An R package with illustrative examples for weighted ATT methods (a generalization to the conventional ATT estimand). The current verison of this package can implment the following estimators: 

* Propensity score (PS) weighting estimator for WATT, with (i) bootstrap; (ii) sandwich variance estimations. For (i), we allow users to specify different PS and outcome regression (OR) models used in `SuperLearner` R package (see <a href="https://academic.oup.com/ectj/article/21/1/C1/5056401">Chernozhukov et al. (2018)</a>) and detailed explanations in the example below. For (ii), we only allow linear model for the OR models and logistic regression for the PS model. Usually, the sandwich variance estimation is designed for parametric models. When using machine learning for PS and OR models, one can also consider a direct variance estimation by the mean square of influence functions, as the bootstrap can be time-consuming. However, for weighting estimator, its influence function under different models is still not clear in literature, so we have not yet included this variance option in the package. 

* Augmented estimator for WATT with only sandwich variance estimator, allowing only linear models for OR and logistic regression for PS. 

## Installation
To install the latest version of the R package from GitHub, please run following commands in R:

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("yiliu1998/ATTweights")
```

In addition, please install and load the following required packages needed for some functions in our package:

```r
if (!require("SuperLearner")) install.packages("SuperLearner")
library(SuperLearner)
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
```

## Use

to be updated...

## Contact 
The R code is maintained by Yi Liu (Please feel free to reach out at  yi.liu.biostat@gmail.com, if you have any questions). 

## Reference
Please cite the following paper:

Liu, Y., Li, H., Zhou, Y., & Matsouaka, R. (2024). Average treatment effect on the treated, under lack of positivity. Statistical Methods in Medical Research. 
