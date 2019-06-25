# take smaller locations
nlocs = 1e3
ncar = 8
nstation = 15
T=1

# estimates from data
load("data_blockmed.Rda")

df_block= df_block_tseg15sec;
# -------
df_block[,colnames(df_block)%in%paste0("PC",1:7)] = scale(df_block[,colnames(df_block)%in%paste0("PC",1:7)])
df_block = df_block[order(df_block$t),]

fooname = c()
foo = paste0("(",paste0("PC",1:7,collapse ="+"),")")
for(i in 1:4) fooname = c(fooname,paste0("H",i,"*",foo))
linmod = lm(formula(paste("Y_block_med ~ ",paste0(fooname,collapse ="+"))),data=df_block,na.action = na.exclude)
df_block$resid = residuals(linmod) #rsq 0.1579

set.seed(1)
require(FNN)
require(fields)
source('Rfunctions_GC.R')

# Select tracks ---------------------------------------------------------
######## set up covariate data ########
covs_name <- "Data_Covariates.csv"
covs      <- read.csv(covs_name,nrows=Inf)

X_names   <- colnames(covs)
junk1     <- grepl("_100",X_names)
junk2     <- grepl("_250",X_names)
junk3     <- grepl("_500",X_names)
junk4     <- grepl("_1000",X_names)
junk5     <- grepl("_2500",X_names)
junk6     <- grepl("_5000",X_names)
junk7     <- grepl("_10K",X_names)
junk      <- junk1 | junk2 | junk3 | junk4 | junk5 | junk6 | junk7
X_sim         <- covs[,!junk]
X_sim <- X_sim[,-c(1:3,6,18:19,21,24:27,39,40)]

pc <- prcomp(X_sim, scale.=TRUE, center=TRUE)
percvar <- cumsum(pc$sdev^2)/sum(pc$sdev^2)

pc_keep <- pc$rotation[,1:7]
X_sim2 <- scale(X_sim)%*%pc_keep

savedatfull <- NULL
datnoon = subset(df_block,hour >= 13 & hour < 16)
keep = tapply(datnoon$hour,floor(datnoon$t),function(x) max(x)==15 & min(x) ==13)
datnoon = subset(datnoon, floor(datnoon$t) %in% names(keep)[keep])
tab  = table(floor(datnoon$t))
tab = tab[(tab>=50)]

for(rep in 1:30){
  randomdays = as.numeric(sample(names(tab),15,replace = F))
  pathdats = subset(datnoon, floor(t) %in% randomdays)
  tapply(pathdats$locID,floor(pathdats$t),length)
  pathlocst = pathdats[,c("locID","Longitude","Latitude","t","hour",paste0("PC",1:7))]
  pathlocst$day = floor(pathdats$t)
  tapply(pathlocst$locID,pathlocst$day,length)
  pathlocst_list <- split(data.frame(pathlocst),pathlocst$day)

  dats = NULL
  day = 1
  for(i in 1:length(pathlocst_list)){
    foo = pathlocst_list[[i]]
    foo$t = foo$t - floor(foo$t) + day
    foo$day = day
    dats = rbind(dats,data.frame(foo, car = paste0("Car",i)))
  }

  dats = dats[order(dats$t),]

  # create station data -----------------------------------------------------
  X_sim2 = data.frame(X_sim2,locID = covs[,1])

  stations = covs[seq(1,nrow(covs),len=nstation),1:3]
  t = seq(12/24,16/24,by=15/24/3600)
  dats_stations = expand.grid(locID=stations[,1], t=t,day = 1)
  dats_stations$hour = floor(dats_stations$t*24)
  dats_stations$t = dats_stations$t+dats_stations$day
  dats_stations$car = "Station"
  dats_stations = merge(dats_stations,stations,by.x ="locID",by.y ="ID")
  dats_stations = merge(dats_stations,X_sim2,by ="locID")
  dats_stations$Longitude = dats_stations$Long30m; dats_stations$Latitude = dats_stations$Lat30m
  dats_stations = dats_stations[order(dats_stations$t),colnames(dats)]

  # create grid for hour 16
  dats_grid = data.frame(locID=covs[,1],Longitude=covs[,2],Latitude=covs[,3], hour=16,day = 1)
  dats_grid$t = dats_grid$hour/24 + dats_grid$day
  dats_grid$car = "grid"
  dats_grid = merge(dats_grid,X_sim2,by ="locID")
  dats_grid = dats_grid[,colnames(dats)]

  # dat full
  datfull = rbind(dats,dats_stations,dats_grid)

  H1         <- sin(2*pi*datfull$hour/24)
  H2         <- cos(2*pi*datfull$hour/24)
  H3         <- sin(4*pi*datfull$hour/24)
  H4         <- cos(4*pi*datfull$hour/24)
  datfull    <- cbind(datfull,H1, H2, H3, H4)
  savedatfull[[rep]] <- datfull
}
save(savedatfull,file="design_deployment.Rda")

# param est ---------------------------------------------------------------
# param estimate from data
# load(file="design_deployment.Rda")
nstation = 15
ncar = 15
params <- c(0.024, 0.043, 0.239, 0.351)
params <- c(0.5492, 0.0188,  0.0049, 0.0038) # 2km and 30 minutes
params <- c(0.3421053*0.452, 0.01918182, 0.1875, (1-0.3421053)*0.452)

# get MSPE
for(rep in 1:length(savedatfull)){
  datfull  = savedatfull[[rep]]
  predcar <- c()
  predst  <- c()
  totalst <- nrow(subset(datfull,car=="grid"))
  dats_grid <- subset(datfull,car=="grid")[seq(1,totalst,l=2000),]
  for(numcar in 1:ncar){
    datcond <- subset(datfull,car%in%paste0("Car",1:numcar) & hour %in% 13:15 & t <= 1+15.5/24)
    predcar[numcar] <- mspe(datcond,dats_grid,
                            as.matrix(datcond[,c("Longitude","Latitude","t")]),
                            as.matrix(dats_grid[,c("Longitude","Latitude","t")]),
                            parms=params)
  }
  
  dats_st = subset(datfull,car=="Station")
  locIDlist = unique(dats_st$locID)
  for(numst in 1:nstation){
    datcond <- subset(dats_st,locID %in% locIDlist[1:numst] & hour %in% 13:15 & t <= 1+15.5/24)
    predst[numst] <- mspe(datcond,dats_grid,
                          as.matrix(datcond[,c("Longitude","Latitude","t")]),
                          as.matrix(dats_grid[,c("Longitude","Latitude","t")]),
                          parms=params)
  }
  
  
  interpcar <- c()
  interpst  <- c()
  dats_grid$t <- 1 + 15.25/24 ## time 14:15pm
  for(numcar in 1:ncar){
    datcond <- subset(datfull,car%in%paste0("Car",1:numcar) & hour %in% 13:15 & t <= 1+15.5/24)
    interpcar[numcar] <- mse(datcond,dats_grid,
                             as.matrix(datcond[,c("Longitude","Latitude","t")]),
                             as.matrix(dats_grid[,c("Longitude","Latitude","t")]),
                             parms=params)
  }
  
  dats_st = subset(datfull,car=="Station")
  locIDlist = unique(dats_st$locID)
  for(numst in 1:nstation){
    datcond <- subset(dats_st,locID %in% locIDlist[1:numst] & hour %in% 13:15 & t <= 1+15.5/24)
    interpst[numst] <- mse(datcond,dats_grid,
                           as.matrix(datcond[,c("Longitude","Latitude","t")]),
                           as.matrix(dats_grid[,c("Longitude","Latitude","t")]),
                           parms=params)
  }
  save(predcar,predst,interpcar,interpst,file=paste0("design_output_composite_rep",rep,".Rda"))
}

# get replicate results
predcarall <- predstall <- interpcarall <- interpstall <- NULL

for(rep in 1:30){
  load(file=paste0("design_output_composite_rep",rep,".Rda"))
  predcarall   <- rbind(predcarall,predcar)
  predstall    <- rbind(predstall,predst)
  interpcarall <- rbind(interpcarall,interpcar)
  interpstall  <- rbind(interpstall,interpst)
}
interpcarall = interpcarall+ (1-0.3421053)*0.452
interpstall = interpstall + (1-0.3421053)*0.452

save(predcarall, predstall, interpcarall, interpstall, file="design_output_all.Rda")
