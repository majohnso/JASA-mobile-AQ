# rm(list=ls())

# require(RANN)
# require(chron)
# require(fields)

# load data ---------------------------------------------------------------
covs_name <- "Data_Covariates.csv" #GIS covariates
goog_name <- "oakland_data_simulated.csv" # raw google car data

covs      <- suppressWarnings(read.csv(covs_name,nrows=Inf))
goog      <- suppressWarnings(read.csv(goog_name,nrows=Inf,as.is=TRUE))

rm(covs_name,goog_name)

goog      <- goog[!is.na(goog[,6]),] # work with NO2 data, rm NAs
goog      <- goog[order(goog$Date_Time),]
s         <- goog[,4:3]
Y         <- goog[,6]
nn        <- as.vector(RANN::nn2(covs[,2:3],s,k=1)$nn.idx)
X         <- as.matrix(covs[nn,-c(1:3)])
Y         <- log(Y)

X_names   <- colnames(X)

## principal component analysis
X_sim         <- covs[,-c(1:3)]
X_sim <- X_sim[,-c(3,15:16,18,21:24,36,37)]

pc <- prcomp(X_sim, scale.=TRUE, center=TRUE)
percvar <- cumsum(pc$sdev^2)/sum(pc$sdev^2)
id <- sum(percvar <= 0.5) + 1
pc_keep <- pc$rotation[,1:7]
X_sim2 <- scale(X_sim)%*%pc_keep
X_PC   <- as.matrix(X_sim2[nn,])

## remove remaining NAs
spID <- covs$ID[nn]
junk      <- is.na(rowSums(cbind(Y,X,s)))
Y         <- Y[!junk]
X         <- X[!junk,]
spID      <- spID[!junk]
X_PC      <- X_PC[!junk,]
goog      <- goog[!junk,]
# start <- lubridate::ymd_hms("2015-01-01 00:00:00",tz="GMT")

pb.txt <- goog[,1]
pb.date <- as.POSIXct(pb.txt, origin = "2015-01-01 00:00:00",tz="America/Los_Angeles")
attributes(pb.date)$tzone <- "America/Los_Angeles"
t <- pb.date

car       <- goog[,2]
s         <- goog[,4:3]
speed     <- goog[,5]
year      <- as.numeric(substr(t,1,4))
month     <- as.numeric(substr(t,6,7))
day       <- as.numeric(substr(t,9,10))
hour      <- as.numeric(substr(t,12,13)) # hour of the day 17 18 19 20 21 22 16 15 23  0  1 14  2
min       <- as.numeric(substr(t,15,16))
sec       <- as.numeric(substr(t,18,19))
wday      <- weekdays(t,abbreviate = T)
date      <- julian(month, day, year,origin=c(1,1,2015))
t         <- julian(month, day, year,origin=c(1,1,2015))+hour/24+min/(60*24)+sec/(60*60*24)
n         <- length(Y)
rm(covs,nn,X_sim,X_sim2)
df <- data.frame(locID = spID, s, speed, Y, Car = car, year, month, day, hour, min, sec, wday, t, X,X_PC)

junk = c("Res_Roads", "NLCD_Barren_50",  
"NLCD_Deciduous_50",  "NLCD_Evergreen_50","NLCD_MixForest_50","NLCD_Shrub_50","NLCD_Herbaceous_50",
"NLCD_Pasture_50","NLCD_Crops_50","NLCD_WoodyWet_50","NLCD_EmergWet_50","Dist2Airport")
df <- df[,!colnames(df)%in%junk]

## add diurnal cycle covariates
df$H1         <- sin(2*pi*df$hour/24)
df$H2         <- cos(2*pi*df$hour/24)
df$H3         <- sin(4*pi*df$hour/24)
df$H4         <- cos(4*pi*df$hour/24)

# save(covs,df,file="df_all.Rda")

######################################################################
#########################  Temporal Aggregation  #####################
######################################################################

## 1 min aggregation
timeseg <- seq(min(df$t)-1*60/(60*60*24),max(df$t)+1*60/(60*60*24),1*60/(60*60*24))
midtseg <- (timeseg[-1] + timeseg[-length(timeseg)])/2
findtseg <- cut(df$t, breaks = timeseg)
levels(findtseg) <- midtseg
findtseg <- as.numeric(as.character(findtseg))
df$tseg1min <- findtseg

## 15sec aggregation
timeseg <- seq(min(df$t)-15/(60*60*24),max(df$t)+15/(60*60*24),15/(60*60*24))
midtseg <- (timeseg[-1] + timeseg[-length(timeseg)])/2
findtseg <- cut(df$t, breaks = timeseg)
levels(findtseg) <- midtseg
findtseg <- as.numeric(as.character(findtseg))
df$tseg15sec <- findtseg

blockmed = function(x,name){
  n = nrow(x)
  if(n%%2==0){idx = which.min(abs(rank(x[,name])-n/2))[1]}else{idx = which.min(abs(rank(x[,name])-(n+1)/2))[1]}
  med = quantile(x[,name],prob=0.5,na.rm=T)
  X = x[idx,]
  return(c(med,X))
}

timelist = c("tseg15sec","tseg1min")
for(tseg in rev(timelist)){
  blockA <- sapply(unique(df[,tseg]),function(x) blockmed(df[which(df[,tseg]==x&df[,"Car"]=="Car_A"),],name="Y"))
  blockB <- sapply(unique(df[,tseg]),function(x) blockmed(df[which(df[,tseg]==x&df[,"Car"]=="Car_B"),],name="Y"))
  fooA <- fooB <- matrix(NA, nrow=ncol(blockA),ncol=nrow(blockA))
  for(i in 1:nrow(blockA)) fooA[,i] <- unlist(blockA[i,])
  for(i in 1:nrow(blockB)) fooB[,i] <- unlist(blockB[i,])
  foo <- rbind(fooA, fooB)
  foo1 <- apply(!is.na(foo),1,all)
  foo <- foo[foo1,]
  colnames(foo) <- rownames(blockA)
  assign(paste0("df_block_",tseg),foo)
  print(paste("Finished:", tseg));
  rm(foo);rm(blockA);rm(blockB)
}

## make it data frame
df_block_tseg15sec <- as.data.frame(df_block_tseg15sec)
df_block_tseg1min <- as.data.frame(df_block_tseg1min)

colnames(df_block_tseg15sec)[1] <- "Y_block_med"
colnames(df_block_tseg1min)[1] <- "Y_block_med"

df_block_tseg15sec <- (df_block_tseg15sec)[,!colnames(df_block_tseg15sec)%in%c("Y")]
df_block_tseg1min <- (df_block_tseg1min)[,!colnames(df_block_tseg1min)%in%c("Y")]

df <- df[order(df$t),]
df_block_tseg15sec <- df_block_tseg15sec[order(df_block_tseg15sec$t),]
df_block_tseg1min <- df_block_tseg1min[order(df_block_tseg1min$t),]
save(df,df_block_tseg15sec,df_block_tseg1min,file="data_blockmed.Rda") # save out to .rda file

