# ATTweights
An R package with illustrative example for weighted ATT methods (a generalization to the conventional ATT). 

## Installation
To install the latest version of the R package from GitHub, please run following commands in R:

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("yiliu1998/ATTweights")
```

In addition, please install and load the following required packages for running function(s) in our package:

```r
if (!require("SuperLearner")) install.packages("SuperLearner")
library(SuperLearner)
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
```

## Demonstration
You can download and run this Rmd file ([click here](https://github.com/yiliu1998/ATTweights/tree/main/vignettes)) in your R Studio after downloading the package, which gives an illustrative example of our package.  

## Contact 
The R code is maintained by Yi Liu (Please feel free to reach out at  yi.liu.biostat@gmail.com, if you have any questions). 

## Reference
Please cite the following paper:

Liu, Y., Li, H., Zhou, Y., & Matsouaka, R. (2024). Average treatment effect on the treated, under lack of positivity. Statistical Methods in Medical Research, accepted. 
