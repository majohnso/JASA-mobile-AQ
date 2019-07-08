README
================

Data Files
----------

### 

R/C++ Files
-----------

### 

Reproducing the analysis on simulated data
------------------------------------------

The R script run\_all\_code.R illustrates the analysis performed in "title" using simulated data. An interested user can obtain the actual data used from [here](https://docs.google.com/forms/d/e/1FAIpQLSf_4GIkK1tmVMFRSxz42KgvOM3Z3NGeOFFje_FS8FBbz1vTig/viewform), in which case the R script should recreate the analysis included in the paper.

The following section walks through the various components of "run\_all\_code.R."

### Set Parameters

``` r
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

### Load R Packages

``` r
pkg.names <- c("Rcpp", "RcppArmadillo", "doParallel", "RANN", "chron", "fields", "dplyr", "FNN", "ggplot2", "ggmap")
for(i in 1:length(pkg.names)){
  if(pkg.names[i] %in% rownames(installed.packages()) == FALSE) {
    install.packages(pkg.names[i])
    cat("\r", pkg.names[i])
  }
}

lp <- sapply(pkg.names, require, character.only = TRUE)

if(sum(lp)!=length(pkg.names)){
  stop(paste0("Unable to load packages: ", paste(pkg.names[!lp], collapse = ", ")))
}
```

### Source R Functions

``` r
source("Rfunctions_GC.R")
sourceCpp("cpp_code_GC.cpp")
```

### Process Raw Data and Perform Temporal Aggregation

``` r
source("data_setup.R")
## Necessary files loaded in data_setup.R:
# Data_IDs.csv # road segment ID
# Data_Covariates.csv #GIS covariates
# oakland_data_simulated.csv # simulated Google car dataset
```

### Model Fitting Setup

``` r
########### Set Up Data For Model Fitting ##########
h <- c(0.02, 5, 15, 60) # minutes
m <- c(10, 30, 60) # minutes 


load("data_blockmed.Rda")
df$Y_block_med <- df$Y
df$Car = ifelse(df$Car=="Car_B",2,1)

covar = c("Longitude","Latitude","t",paste0("PC",1:7))

if(type == 1) df_block = df
if(type == 2) df_block = df_block_tseg15sec
if(type == 3) df_block = df_block_tseg1min

keep = c("Y_block_med","locID","Longitude","Latitude","Car","speed","year",
         "month","day","hour","min","sec","wday","t",paste0("PC",1:7),paste0("H",1:4))
df_block[,colnames(df_block)%in%paste0("PC",1:7)] = scale(df_block[,colnames(df_block)%in%paste0("PC",1:7)])
df_block=df_block[order(df_block$t),keep]

regressorname = c()
foo = paste0("(",paste0("PC",1:7,collapse ="+"),")")
for(i in 1:4) regressorname = c(regressorname,paste0("H",i,"*",foo))

PCidx <- which(colnames(df_block)%in%c(paste0("PC",1:7),paste0("H",1:4)))
cov_names <- c(paste0("PC",1:7),paste0("H",1:4))

days = sort(unique(floor(df_block$t)))
df_block$day = floor(df_block$t)

if(rolling_window){
  #  create blocks of two weeks
  wk_block = seq(192,max(days),by = 7*1)
  daysint = findInterval(df_block$day, wk_block)
  df_block$daysint = daysint; 
  df_block$idx = 1:nrow(df_block);
}
```

### Fit Models

``` r
source("st_stx_estimation_prediction.R") 
# ST and STx Vecchia estimation and prediction 
# for type, and combination of h (j) and m (jj)
# NOTE: this takes a long time to run

if(rolling_window){
  source("rolling_window_estimation.R")
}

source("spatial_only_estimation_prediction.R")
# S model Vecchia approximation estimation and prediction
# for type, and combination of h (j) and m (jj)
# NOTE: this takes a long time to run
```

### Map Forecasts

``` r
# run for dayidx = 253 and dayidx = 490 
if(type == 2){
  for(dayidx in c(253, 490)){
    source("15min_map_forecasts.R")
  }
}
```

### Mobile vs. Stationary Simulation

``` r
source("deployment_design.R")
# NOTE: this takes a long time to run
```