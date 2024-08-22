# Research Code
Implementation code for the article on non-parametric estimation of transition intensities in interval censored Markov multi-state models without loops.


This repository aims to make the results from [Non-parametric estimation of transition intensities in inverval-censored Markov multi-state models without loops](www.google.com) reproducible. 
There are two directories, with `Data Analysis` showing the code used to obtain results for the [Signal-Tandmobiel](https://rdrr.io/cran/icensBKL/man/tandmobAll.html) data and `Simulations` showing the code for the simulation study.

The code in this repository makes use of the `icmstate` package, which can be installed as follows:

```{r}
library(devtools)
install_github("d-gomon/icmstate")
```


The two subdirectories are described in more detail below.

## Data Analysis

In this directory the code is contained which was used to perform analyses on the [Signal-Tandmobiel](https://rdrr.io/cran/icensBKL/man/tandmobAll.html) study. The results can be reproduced exactly, and should not take too long to run on an average consumer PC. Click the directory for more information.


## Simulations

In this directory the code used to perform the simulation studies can be found. Although exactly reproducible, the simulation study was run externally on a computing cluster, so attempting to run the code on a consumer PC is not recommended. Click the directory for more information.

