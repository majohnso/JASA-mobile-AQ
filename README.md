---
title: "README"
date: "7/8/2019"
output: html_document
---

## Data Files

###

## R/C++ Files
###

## Reproducing the analysis on simulated data

The R script run_all_code.R illustrates the analysis performed in "title" using simulated data. An interested user can obtain the actual data used from 'website', in which case the R script should recreate the analysis included in the paper. 

The following section walks through the various components of "run_all_code.R."

### Analysis Parameters
```{r setup_params, eval=FALSE}

rolling_window = FALSE # should rolling window estimation be performed (computationally expensive)

type <- 2; # type in 1:3, defines which data product to use 
# 1 == raw data, 2 == 15sec aggregates, 3 == 1min aggregates

j <- 3 # j in 1:4, which estimation lag, index of h <- c(0.02, 5, 15, 60) min
jj <- 3 # jj in 1:3, which size of conditioning set, index of m <- c(10, 30, 60) min

lag <- 21 # lag for moving window 

if(type==1){
  ns <- 100 # if using 1sec data, define the size to subsample the conditioning set in Vecchia
  nsA <- 800 # if using 1sec data, define the nearest neighbor size for Car A prediction
}

nnbs <- 800 # number of nearest neighbors for spatial only prediction

maxit <- 1e3 # maximum number of iterations for Nelder-Mead maximization of Vecchia log-likelihood
```


