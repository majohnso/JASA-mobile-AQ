README
================

This repository contains simulated data, covariates, and code to illustrate the methodology and analysis published in \`\`Fine-scale spatiotemporal air pollution analysis using mobile monitors on Google Street View vehicles" (add DOI).

Data Files
----------

#### oakland\_data\_simulated.csv

-   Contains simulated Oakland Air Quality data
-   Variables include:
    -   <tt>Date\_Time</tt>: Data and time (seconds) of each measurement
    -   <tt>Car\_Identifier</tt>: identifies which car (A or B) obtained the measurement
    -   <tt>Latitude</tt>: latitude of car position at time of measurement
    -   <tt>Longitude</tt>: longitude of car position at time of measurement
    -   <tt>Car\_Speed</tt>: speed of car (km/s)
    -   <tt>NO2</tt>: measured nitrogen dioxide (ppm)

##### Simulation Details

Data are simulated for a subset of 60 days between Sept. 11, 2015 and May 5, 2016. Spatial locations and measurement times are sampled from actual observed car paths, with hour of observation altered to build a sequence of data from 9:00 to 18:00 on each day. Measurements of N02 at the simulated paths and times are constructed as white noise perturbed predictions from the fitted ST model estimated using the actual Google air quality data. The full sequence of simulated data is then subsampled to approximately 1/6 the original size to produce 76,484 total simulated measurements, roughly 10% the size of the actual Oakland data.

Note: the simulated data is meant to help understand how the following code can be used but is not intended to show any interesting results. A interested user should obtain the actual data from Google. The data is freely available upon request [here](https://docs.google.com/forms/d/e/1FAIpQLSf_4GIkK1tmVMFRSxz42KgvOM3Z3NGeOFFje_FS8FBbz1vTig/viewform).

#### Data\_Covariates.csv

-   Contains (real) spatial covariate data for 30m road segments in Oakland.
-   Variables include:
    -   <tt>ID</tt>: unique identifier for each 30m road segment
    -   <tt>Long30m</tt>: longitude at the center of the segment
    -   <tt>Lat30m</tt>: latitude at the center of the segment
    -   <tt>Elevation\_50</tt>: mean elevation within a circular buffer of 50 meters
    -   <tt>Total\_Road\_50</tt>: total length of highways, major road, residential roads within a circular 50 meter buffer region
    -   <tt>Hwy\_Road\_50</tt>: highways road length within a circular buffer of 50 meters
    -   <tt>Maj\_Road\_50</tt>: major arterials road length within a circular buffer of 50 meters
    -   <tt>Res\_Road\_50</tt>: residential roads within a circular buffer of 50 meters
    -   <tt>Pop\_50</tt>: population density (people/sq-km) within a circular buffer 50 meters
    -   Binary road classification variables:
        -   <tt>Hwy\_Roads</tt>: highways
        -   <tt>Major\_Roads</tt>: major through road
        -   <tt>Res\_Roads</tt>: residential road
        -   <tt>Local\_Trucks</tt>: designated heavy-duty truck routes
        -   <tt>Local\_Restricted\_Trucks</tt>: restricted heavy-duty truck routes
        -   <tt>Commerical.Zone</tt>: commercial zoning
        -   <tt>Industrial.Zone</tt>: industrial zoning
        -   <tt>Residential.Zone</tt>: residential zoning
        -   <tt>NDVI\_50</tt>: The average Normalized Difference Vegetative Index within circular buffer of 50 meters
    -   Variables constructed based on the National Land Cover Database satellite imagery file. Each variable shows the percent of land cover of type within a circular buffer of 50 meters.
        -   <tt>NLCD\_Water\_50</tt>: water
        -   <tt>NLCD\_DevOpen\_50</tt>: open areas in developments
        -   <tt>NLCD\_DevLow\_50</tt>: low developed areas
        -   <tt>NLCD\_DevMed\_50</tt>: medium developed areas
        -   <tt>NLCD\_DevHigh\_50</tt>: highly developed areas
        -   <tt>NLCD\_Barren\_50</tt>: barren
        -   <tt>NLCD\_Deciduous\_50</tt>: deciduous forest
        -   <tt>NLCD\_Evergreen\_50</tt>: evergreen forest
        -   <tt>NLCD\_MixForest\_50</tt>: mixed forest
        -   <tt>NLCD\_Shrub\_50</tt>: shrublands
        -   <tt>NLCD\_Herbaceous\_50</tt>: herbaceous
        -   <tt>NLCD\_Pasture\_50</tt>: pasture
        -   <tt>NLCD\_Crops\_50</tt>: crops
        -   <tt>NLCD\_WoodyWet\_50</tt>: woody wetlands
        -   <tt>NLCD\_EmergWet\_50</tt>: emerging wetlands
        -   <tt>NLCD\_Impervious\_50</tt>: impervious surface
    -   Variables representing cumulative exponentially decaying contribution from point sources
        -   <tt>Distance\_to\_NPL</tt>: National Priority Listing (NPL) sites
        -   <tt>Distance\_to\_Rail</tt>: railroads
        -   <tt>Distance\_to\_TRI</tt>: Toxic Release Inventory (TRI) sites
    -   Variables representing minimum inverse distance to point sources
        -   <tt>MinDist2Port</tt>: all principal port and facility locations
        -   <tt>MinDist2MainPort</tt>: principal port locations
        -   <tt>Dist2MainAirport</tt>: Major airports
        -   <tt>Dist2Airport</tt>: airports

R/C++ Files
-----------

<tt> Rfunctions\_GC.R </tt>

-   contains all R functions to perform estimation and prediction for the STx, ST and S models

<tt> cpp\_code\_GC.cpp </tt>

-   contains helper c++ files to speed up estimation and prediction

<tt> run\_all\_code.R </tt>

-   <tt> R </tt> script to perform analysis included in the paper on simulated data provided in the Data folder
-   Calls the follwing helper scripts
    -   <tt> data\_setup.R </tt>
        -   cleans up raw data and performs temporal aggregation at 15 sec and 1 min block medians
        -   each aggregated dataset is saved and stored in <tt> data\_blockmed.Rda </tt>
    -   <tt> st\_stx\_estimation\_prediction.R </tt>
        -   performs parameter estimation and prediction for the ST and STx models
        -   results are stored as <tt> R </tt> objects and saved out in <tt> .Rda </tt> files
    -   <tt> spatial\_only\_estimation\_prediciton.R </tt>
        -   performs parameter estimation and prediction for the S (spatial only) model
        -   results are stored as <tt> R </tt> objects and saved out in <tt> .Rda </tt> files
    -   <tt> 15min\_map\_forecasts.R </tt>
        -   creates 15 min ahead spatial map forecasts for the ST model for two days
        -   used to recreate maps contained in Figures 5 and 6 of the paper
    -   <tt> rolling\_window\_estimation.R </tt>
        -   performs ST and STx model estimation on data in 21-week rolling windows
        -   results are stored in <tt> .Rda </tt> files
    -   <tt> deployment\_design.R </tt>
        -   performs the simulation experiment comparing mspe for mobile vs fixed-location monitors

<tt> plotting\_code.R </tt>

-   Contains code to construct the plots include in the paper after obtaining all results from <tt>run\_all\_code.R</tt>

Reproducing the analysis on simulated data
------------------------------------------

The R script <tt> run\_all\_code.R </tt> illustrates the analysis performed in \`\`Fine-scale spatiotemporal air pollution analysis using mobile monitors on Google Street View vehicles" using simulated data. If the user obtains the actual Google data, <tt> run\_all\_code.R </tt> should recreate the analysis included in the paper.

The following sections walk through the various components of <tt> run\_all\_code.R </tt>. Note that this <tt>R</tt> script and all others assume the data and scripts are all stored in a single, common folder.

### Set Parameters

``` r
rolling_window = FALSE # should rolling window estimation be performed? (computationally expensive)

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
pkg.names <- c("Rcpp", "RcppArmadillo", "doParallel", "RANN", "chron", "fields", 
               "dplyr", "FNN", "ggplot2", "ggmap")

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
#
# Saves out temporally aggregated datasets in "data_blockmed.Rda"
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
# NOTE: this estimation can take several hours to run, depending on the number of cores used
```

ST and STx Vecchia estimation and prediction for <tt>type</tt>, and combination of <tt>h</tt> (<tt>j</tt>) and <tt>m</tt> (<tt>jj</tt>) specified above. Saves out estimated parameters; 5min, 15min and 60min forecasts; Car A predictions; mspe and correlation in <tt>.Rda</tt> files.

``` r
source("spatial_only_estimation_prediction.R")
# NOTE: this can take several hours to run, depending on the number of cores used.
```

S model Vecchia approximation estimation and prediction for <tt>type</tt>, and combination of <tt>h</tt> (<tt>j</tt>) and <tt>m</tt> (<tt>jj</tt>) specified above. Saves out estimated parameters, spatial predictions, Car A predictions, mspe and correlation in <tt>.Rda</tt> files

``` r
if(rolling_window){
  source("rolling_window_estimation.R")
}
```

ST and STx Vecchia estimation and prediction using rolling windows with window size lag = 21 weeks. For example, data from week 1-21 are used for parameter estimation, then prediction is made for week 22. This procedure is then repeated for the next window of week 2-22. The rolling window analysis to compare ST and STx Vecchia models is performed in parallel for each window in practice. It is not performed here by default (<tt>rolling\_window = FALSE</tt>) because the example simulated data spans only two months. The results are saved out for each window in <tt>.Rda</tt> files

### Map Forecasts

``` r
# run for dayidx = 253 and dayidx = 490 
if(type == 2){
  for(dayidx in c(253, 490)){
    source("15min_map_forecasts.R")
  }
}
```

Creates full spatial 15-min ahead map forecasts for the two different days and times included in the paper (Figures 5 and 6).

### Mobile vs. Stationary Simulation

``` r
source("deployment_design.R")

# NOTE: this can take several hours to run, depending on the number of cores used.
```

Performs a simulation experiment comparing mean squared prediction error (MSPE) for mobile vs fixed-location monitors for short-term forecasting and spatial interpolation. For <tt>ncar=1,...,15</tt> and <tt>nstation=1,...,15</tt>, we sample <tt>ncar</tt> number of routes from the Google data and randomly select <tt>nstation</tt> number of locations in the study region as fixed-location monitors. The MSPE is computed for both type of monitors conditional on locations from ncar and nstation of monitors respectively. This procedure is repeated 30 times to obtain uncertainty of the MSPE.
