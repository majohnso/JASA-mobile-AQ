##########################################################################
#################     Static Estimation & Prediction    ##################
##########################################################################

#### build training data
est_data <- subset(df_block,floor(t) %in% 301:351) #P2 total var = 0.4685434

# est_data <- subset(df_block,floor(t) < 351)
byday <- split(est_data, floor(est_data$t))
only_after_2hrs <- lapply(byday, function(u) subset(u, t > min(t) + 2/24))

est_data2 <- do.call(rbind, only_after_2hrs)
idx <- which(est_data$t %in% est_data2$t)

locstime <- est_data[,c("Longitude","Latitude","t")]
locstimeX <- est_data[,c("Longitude","Latitude","t",paste0("PC",1:7))]

linmod = lm(formula(paste("Y_block_med ~ ",paste0(regressorname,collapse ="+"))),data=est_data,na.action = na.exclude)
est_data$resid = residuals(linmod) 
y <- est_data$resid
total_variance = var(residuals(linmod))
save(linmod, file=paste0("linmod_type", type, ".Rda"))

################### Vecchia for space-time and space-time-x parameters #################
### space-time model
startparms <- c(9,0.01,0.01) 
if(type!=1){
  fit_st <- optim(log(startparms), fun_to_max_vecchia.st, var = total_variance,y=y, idx=idx, locstime=locstime, h=h[j], m=m[jj], ncore=ncore, control=list(maxit=maxit, trace=10))
}else{
  fit_st <- optim(log(startparms), fun_to_max_vecchia.st_smple, var = total_variance,y=y, idx=idx, locstime=locstime, h=h[j], m=m[jj], ncore=ncore, nsize=ns, control=list(maxit=maxit, trace=10))
}

### space-time-x model
startparms <- c(9,0.01,0.01,2) 
if(type!=1){
  fit_stx <- optim(log(startparms), fun_to_max_vecchia.stx, var = total_variance,y=y, idx=idx, locstime=locstimeX, h=h[j], m=m[jj], ncore=ncore, control=list(maxit=maxit, trace=10))
}else{
  fit_stx <- optim(log(startparms), fun_to_max_vecchia.stx_smple, var = total_variance,y=y, idx=idx, locstime=locstimeX, h=h[j], m=m[jj], ncore=ncore, nsize=ns, control=list(maxit=maxit, trace=10))
}

cat("\n")

################### Prediction ###############
datpred <- subset(df_block,floor(t) %in% 352:400) #P2
byday <- split(datpred, floor(datpred$t))
only_after_1hrs <- lapply(byday, function(u) subset(u, t > min(t) + 1/24))

datpred2 <- do.call(rbind, only_after_1hrs)
id_pred <- which(datpred$t %in% datpred2$t)

parms_st = exp(fit_st$par)
parms_st[1] = parms_st[1]/(parms_st[1]+1)
parms_st = c(parms_st[1]*total_variance,parms_st[2],parms_st[3],1,(1-parms_st[1])*total_variance)

parms_stx = exp(fit_stx$par)
parms_stx[1] = parms_stx[1]/(parms_stx[1]+1)
parms_stx = c(parms_stx[1]*total_variance,parms_stx[2:4],(1-parms_stx[1])*total_variance)

for(hh in c(5,15,60)){
  if(type!=1){
    dat1 <- mclapply(id_pred, function(x) krig.mh.st(x, datpred, parms_st, h=hh, m=60, cov_names, linmod), mc.cores = ncore)
  }else{
    dat1 <- mclapply(id_pred, function(x) krig.mh.st_smple(x, datpred, parms_st, h=hh, m=60, cov_names, linmod, nsize=ns), mc.cores = ncore)
  }
  dat1 <- as.data.frame(do.call(rbind, dat1))
  names(dat1) <- c("Pred_St", "Var_St", "Pred_Xb")
  
  if(type!=1){
    dat2 <- mclapply(id_pred, function(x) krig.mh.stx(x, datpred, parms_stx, h=hh, m=60, cov_names, linmod), mc.cores = ncore)
  }else{
    dat2 <- mclapply(id_pred, function(x) krig.mh.stx_smple(x, datpred, parms_stx, h=hh, m=60, cov_names, linmod, nsize=ns), mc.cores = ncore)
  }
  dat2 <- as.data.frame(do.call(rbind, dat2))
  names(dat2) <- c("Pred_Stx", "Var_Stx", "Pred_Xb")
  
  dat <- cbind(dat1[,1:2], dat2)
  
  dat$Y <- datpred$Y_block_med[id_pred]
  dat$idx <- datpred$idx[id_pred]
  dat <- dat %>% mutate(
    resid_st = Y - Pred_St,
    resid_stx = Y - Pred_Stx,
    resid_xb = Y - Pred_Xb
  )
  
  mspes <- apply(dat[,c("resid_st","resid_stx","resid_xb")],2,MSE)
  cors <- apply(dat[,c("Pred_St","Pred_Stx","Pred_Xb")],2,function(x) cor(x,dat[,"Y"]))
  assign(paste0("pred_",hh,"min"), dat)
  assign(paste0("mspe_",hh,"min"), mspes)
  assign(paste0("cors_",hh,"min"), cors)
  
  cat("\r",paste0("Finished: ST/STx ",hh," min ahead forecasts\n"))
}

mspes = rbind(mspe_5min, mspe_15min, mspe_60min)
cors =  rbind(cors_5min, cors_15min, cors_60min)

save(linmod,pred_5min,pred_15min,pred_60min,mspes,cors,fit_st,fit_stx,
     file=paste0("composite_static_lambda_type",type,"h",h[j],"_m",m[jj],".Rda"))

############### Car A Prediction ################
datpred <- subset(df_block,floor(t) > 351) #P2
byday <- split(datpred, floor(datpred$t))
kp <- as.numeric(names(which(sapply(byday, function(u) sum(u$Car==1) > 0 & sum(u$Car==2) > 0))))

datpred <- subset(df_block,floor(t) %in% kp) #P2
byday <- split(datpred, floor(datpred$t))
only_after_1hrs <- lapply(byday, function(u) subset(u, hour > min(hour)))
datpred$hr_day <- floor(datpred$t) + datpred$hour/24

datpred2 <- do.call(rbind, only_after_1hrs)
datpred2$hr_day <- floor(datpred2$t) + datpred2$hour/24
hr_pred <- unique(datpred$hr_day[datpred$hr_day %in% datpred2$hr_day])
id_pred <- which(datpred$hr_day %in% datpred2$hr_day)

###### ST prediction
parms= exp(fit_st$par)
parms[1] = parms[1]/(parms[1]+1)
parms = c(parms[1]*total_variance,parms[2:4],(1-parms[1])*total_variance)

if(type!=1){
  dat <- mclapply(hr_pred, function(x){krig.mh.st_carA(x, datpred, parms, cov_names, linmod)}, mc.cores = ncore)
}else{
  dat <- mclapply(hr_pred, function(x){krig.mh.st_carA_smple(x, datpred, parms, cov_names, linmod, nsize=nsA)}, mc.cores = ncore)
}
dat <- as.data.frame(do.call(rbind, dat))
names(dat) <- c("Pred_St", "Var_St", "Pred_Xb")
dat$Y <- subset(datpred[id_pred,], Car==1)$Y_block_med

dat <- dat %>% mutate(
  resid_st = Y - Pred_St,
  resid_xb = Y - Pred_Xb
)

pred_carA <- dat

mspes_carA <- apply(dat[,c("resid_st","resid_xb")],2,MSE)
cors_carA <- apply(dat[,c("Pred_St","Pred_Xb")],2,function(x) cor(x,dat[,"Y"]))

save(mspes_carA, cors_carA, pred_carA,
     file=paste0("CarA_preds_ST_type",type,"h",h[j],"_m",m[jj],".Rda"))

cat("\r",paste0("Finished: ST Car A Predictions\n"))

##### STX car A prediction
parms <- exp(fit_stx$par)
parms[1] = parms[1]/(parms[1]+1)
parms = c(parms[1]*total_variance,parms[2:4],(1-parms[1])*total_variance)

if(type!=1){
  dat <- mclapply(hr_pred, function(x){krig.mh.stx_carA(x, datpred, parms, cov_names, linmod)}, mc.cores = ncore)
}else{
  dat <- mclapply(hr_pred, function(x){krig.mh.stx_carA_smple(x, datpred, parms, cov_names, linmod, nsize=nsA)}, mc.cores = ncore)
}
dat <- as.data.frame(do.call(rbind, dat))
names(dat) <- c("Pred_Stx", "Var_Stx", "Pred_Xb")
dat$Y <- subset(datpred[id_pred,], Car==1)$Y_block_med

dat <- dat %>% mutate(
  resid_stx = Y - Pred_Stx,
  resid_xb = Y - Pred_Xb
)

pred_carA <- dat

mspes_carA <- apply(dat[,c("resid_stx","resid_xb")],2,MSE)
cors_carA <- apply(dat[,c("Pred_Stx","Pred_Xb")],2,function(x) cor(x,dat[,"Y"]))

save(mspes_carA, cors_carA, pred_carA,
     file=paste0("CarA_preds_STX_type",type,"h",h[j],"_m",m[jj],".Rda"))

cat("\r",paste0("Finished: STx Car A Predictions\n"))



