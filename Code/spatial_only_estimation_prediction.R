##########################################################################
#################        Estimation & Prediction        ##################
##########################################################################

################### Vecchia for spatial parameters #################
est_data <- subset(df_block,floor(t) %in% 301:351) #P2 total var = 0.4685434

# est_data <- subset(df_block,floor(t) < 351)
byday <- split(est_data, floor(est_data$t))
only_after_2hrs <- lapply(byday, function(u) subset(u, t > min(t) + 2/24))

est_data2 <- do.call(rbind, only_after_2hrs)
idx <- which(est_data$t %in% est_data2$t)

locstime <- est_data[,c("Longitude","Latitude","t")]

linmod = lm(formula(paste("Y_block_med ~ ",paste0(regressorname,collapse ="+"))),data=est_data,na.action = na.exclude)
est_data$resid = residuals(linmod) 
y <- est_data$resid
total_variance = var(residuals(linmod))

startparms <- c(9,0.01)
if(type != 1){
  fit_s <- optim(log(startparms), fun_to_max_vecchia.s, var = total_variance, y=y, idx=idx, locstime=locstime, m=60, ncore=ncore, control=list(maxit=maxit, trace=10))
}else{
  fit_s <- optim(log(startparms), fun_to_max_vecchia.s_smple, var = total_variance, y=y, idx=idx, locstime=locstime, m=60, ncore=ncore, nsize=ns, control=list(maxit=maxit, trace=10))
}

cat("\n", "Finished: Spatial Parameter Estimation\n")

################### Prediction ###############
datpred <- subset(df_block,floor(t) %in% 352:400) #P2
byday <- split(datpred, floor(datpred$t))
only_after_1hrs <- lapply(byday, function(u) subset(u, t > min(t) + 1/24))

datpred2 <- do.call(rbind, only_after_1hrs)
id_pred <- which(datpred$t %in% datpred2$t)

parms_s = exp(fit_s$par)
parms_s[1] = parms_s[1]/(parms_s[1]+1)
parms_s = c(parms_s[1]*total_variance,parms_s[2],1,1,(1-parms_s[1])*total_variance)

cov_names <- c(paste0("PC",1:7),paste0("H",1:4))
nbs <- nn2(est_data[,c("Longitude","Latitude")], datpred[,c("Longitude","Latitude")],k=800)$nn.idx


dat <- mclapply(id_pred, function(x) krig.mh.s(x, nbs, datpred, est_data, cov_names, parms_s, linmod), mc.cores = ncore)
dat <- as.data.frame(do.call(rbind, dat))
names(dat) <- c("Pred_S","Var_S","Pred_Xb")
pred_s <- dat
pred_s$Y <- datpred[id_pred,]$Y_block_med

mspes <- mean((pred_s$Pred_S - datpred[id_pred,]$Y_block_med)^2)
cors <- cor(pred_s$Pred_S, datpred[id_pred,]$Y_block_med)

if(type==1){
  save(pred_s, mspes, cors, fit_s,
       file=paste0("nbs800_spatial_type",type,"ns",ns,".Rda"))
}else{
  save(pred_s, mspes, cors, fit_s,
       file=paste0("nbs800_spatial_type",type,".Rda"))
}

cat("\r", "Finished: Spatial Forecasts\n")

######################### Car A Prediction ######################
datpred <- subset(df_block,floor(t) > 351) #P2
byday <- split(datpred, floor(datpred$t))
kp <- as.numeric(names(which(sapply(byday, function(u) sum(u$Car==1) > 0 & sum(u$Car==2) > 0))))

datpred <- subset(df_block,floor(t) %in% kp) #
byday <- split(datpred, floor(datpred$t))
only_after_1hrs <- lapply(byday, function(u) subset(u, t > min(t) + 1/24))

datpred2 <- do.call(rbind, only_after_1hrs)
id_pred <- which(datpred$t %in% datpred2$t)
datpred <- subset(datpred[id_pred,], Car==1)

nbs <- nn2(est_data[,c("Longitude","Latitude")], datpred[,c("Longitude","Latitude")],k=nnbs)$nn.idx

dat <- mclapply(1:nrow(datpred), function(x) krig.mh.s(x, nbs, datpred, est_data, cov_names, parms_s, linmod), mc.cores = ncore)
dat <- as.data.frame(do.call(rbind, dat))
names(dat) <- c("Pred_S","Var_S","Pred_Xb")
pred_CarA <- dat
pred_CarA$Y <- datpred$Y_block_med

mspes_CarA <- mean((pred_CarA$Pred_S - datpred$Y_block_med)^2)
cors_CarA <- cor(pred_CarA$Pred_S, datpred$Y_block_med)

save(pred_CarA, mspes_CarA, cors_CarA, 
     file=paste0("nbs",nnbs,"_spatial_CarA_type",type,".Rda"))

cat("\r", "Finished: Spatial Car A Forecasts\n")
