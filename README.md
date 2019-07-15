README
================

Data Files
----------

#### oakland\_data\_simulated.csv

-   Contains simulated Oakland Air Quality data
-   Variables include:
    -   Date\_Time
    -   Car\_Identifier
    -   Latitude
    -   Longitude
    -   Car\_Speed
    -   NO2

#### Data\_IDs.csv

-   Contains 30m road segment IDs
-   Variables include:
    -   Longitude
    -   Latitude
    -   ID

#### Data\_Covariates.csv

-   Contains spatial covariate data
-   Variables include:
    -   ID
    -   Long30m
    -   Lat30m
    -   Hwy\_Roads
    -   Major\_Roads
    -   Res\_Roads
    -   Local\_Trucks
    -   Local\_Restricted\_Trucks
    -   Commerical.Zone
    -   Industrial.Zone
    -   Residential.Zone
    -   NDVI\_50
    -   NLCD\_Water\_50
    -   NLCD\_DevOpen\_50
    -   NLCD\_DevLow\_50
    -   NLCD\_DevMed\_50
    -   NLCD\_DevHigh\_50
    -   NLCD\_Barren\_50
    -   NLCD\_Deciduous\_50
    -   NLCD\_Evergreen\_50
    -   NLCD\_MixForest\_50
    -   NLCD\_Shrub\_50
    -   NLCD\_Herbaceous\_50
    -   NLCD\_Pasture\_50
    -   NLCD\_Crops\_50
    -   NLCD\_WoodyWet\_50
    -   NLCD\_EmergWet\_50
    -   Impervious\_50
    -   Distance\_to\_NPL
    -   Distance\_to\_Rail
    -   Elevation\_50
    -   Distance\_to\_TRI
    -   Total\_Road\_50
    -   Hwy\_Road\_50
    -   Maj\_Road\_50
    -   Res\_Road\_50
    -   Pop\_50
    -   MinDist2Port
    -   MinDist2MainPort
    -   Dist2MainAirport
    -   Dist2Airport

R/C++ Files
-----------

-   <tt> data\_setup.R </tt>
-   <tt> cpp\_code\_GC.cpp </tt>
-   <tt> Rfunctions\_GC.R </tt>
-   <tt> rolling\_window\_estimation.R </tt>
-   <tt> spatial\_only\_estimation\_prediciton.R </tt>
-   <tt> st\_stx\_estimation\_prediction.R </tt>
-   <tt> deployment\_design.R </tt>
-   <tt> 15min\_map\_forecasts.R </tt>
-   <tt> plotting\_code.R </tt>
-   <tt> run\_all\_code.R </tt>

Reproducing the analysis on simulated data
------------------------------------------

The R script <tt> run\_all\_code.R </tt> illustrates the analysis performed in \`\`Fine-scale spatiotemporal air pollution analysis using mobile monitors on Google Street View vehicles" using simulated data. An interested user can obtain the actual data used from [here](https://docs.google.com/forms/d/e/1FAIpQLSf_4GIkK1tmVMFRSxz42KgvOM3Z3NGeOFFje_FS8FBbz1vTig/viewform), in which case the <tt> R </tt> script should recreate the analysis included in the paper.

The following section walks through the various components of <tt> run\_all\_code.R </tt>.

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
# Necessary files loaded in data_setup.R:
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
