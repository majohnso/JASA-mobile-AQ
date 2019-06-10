################################ 15min Map Forecasts #############################
load("data_blockmed.Rda")
df$Y_block_med <- df$Y
df$Car = ifelse(df$Car=="Car_B",2,1)

covar = c("Longitude","Latitude","t",paste0("PC",1:7))

regressorname = c()
foo = paste0("(",paste0("PC",1:7,collapse ="+"),")")
for(i in 1:4) regressorname = c(regressorname,paste0("H",i,"*",foo))

df_block = df
ID = df_block$locID
ID_idx = which(!duplicated(ID))
PCidx <- which(colnames(df_block)%in%c(paste0("PC",1:7),paste0("H",1:4)))

pred_block= df_block[ID_idx,c(1:3,PCidx)]
save(pred_block, file="pred_block.Rda")

nh <- sum(unique(df_block$hour[floor(df_block$t)==dayidx]) > 11)

fore_dat = do.call(rbind, rep(list(pred_block),nh))
fore_dat$hour <- rep(c(12:(12+nh-1)), each=nrow(pred_block))
fore_dat = fore_dat[order(fore_dat$hour),]

fore_dat$H1         <- sin(2*pi*fore_dat$hour/24)
fore_dat$H2         <- cos(2*pi*fore_dat$hour/24)
fore_dat$H3         <- sin(4*pi*fore_dat$hour/24)
fore_dat$H4         <- cos(4*pi*fore_dat$hour/24)
fore_dat$t <- dayidx + fore_dat$hour/24

# type <- 2
df_block = df_block_tseg15sec
df_block[,colnames(df_block)%in%paste0("PC",1:7)] = scale(df_block[,colnames(df_block)%in%paste0("PC",1:7)])
df_block=df_block[order(df_block$t),]

cov_names <- c(paste0("PC",1:7),paste0("H",1:4))

est_data <- subset(df_block,floor(t) %in% 301:351) #P2 total var = 0.4685434

byday <- split(est_data, floor(est_data$t))
only_after_2hrs <- lapply(byday, function(u) subset(u, t > min(t) + 2/24))

est_data2 <- do.call(rbind, only_after_2hrs)
idx <- which(est_data$t %in% est_data2$t)

locstime <- est_data[,c("Longitude","Latitude","t")]
linmod = lm(formula(paste("Y_block_med ~ ",paste0(regressorname,collapse ="+"))),data=est_data,na.action = na.exclude)
est_data$resid = residuals(linmod) 
y <- est_data$resid
total_variance = var(residuals(linmod))

##################### NEW SPATIAL PREDICTION
load(paste0("nbs",nnbs,"_spatial_type",type,".Rda"))

parms_s = exp(fit_s$par)
parms_s[1] = parms_s[1]/(parms_s[1]+1)
parms_s = c(parms_s[1]*total_variance,parms_s[2],1,1,(1-parms_s[1])*total_variance)

sp_data <- fore_dat[,c("Longitude", "Latitude", cov_names)]
nbs <- nn2(est_data[,c("Longitude","Latitude")], sp_data[,1:2],k=240)$nn.idx
sp_dat <- mclapply(1:nrow(sp_data), function(x) krig.mh.s(x, nbs, sp_data, est_data, cov_names, parms_s, linmod), mc.cores = ncore)
sp_dat <- as.data.frame(do.call(rbind, sp_dat))
names(sp_dat) <- c("Pred_S","Var_S","Pred_Xb")
sp_dat$hour <- rep(12:(12+nh-1), each=nrow(pred_block))

save(sp_dat, file=paste0("s_map_preds", dayidx,".Rda"))

##################### space-time prediction
# load("composite_static_lambda_type2h15_m60.Rda")
load(paste0("composite_static_lambda_type",type,"h",h[j],"_m",m[jj],".Rda"))

parms = exp(fit_st$par)
parms[1] = parms[1]/(parms[1]+1)
parms = c(parms[1]*total_variance,parms[2],parms[3],1,(1-parms[1])*total_variance)

results_hour_st <- rep(list(NA),nh)
for(i in 1:nh){
  dat_hour <- subset(fore_dat, hour == c(12:(12+nh-1))[i])
  kp <- which(df_block$t < min(dat_hour$t) & floor(df_block$t) > dayidx - 2)
  
  cond_dat <- df_block[kp,c("Longitude", "Latitude", "t", cov_names, "Y_block_med")]
  datpred <- rbind(cond_dat,cbind(dat_hour[,c("Longitude", "Latitude", "t", cov_names)],Y_block_med=NA))
  id_pred <- (nrow(cond_dat)+1):nrow(datpred)
  
  #### get predictions
  dat_st <- mclapply(id_pred, function(x) krig.mh.st(x, datpred, parms, 15, m=60, cov_names, linmod), mc.cores = ncore)
  dat_st <- as.data.frame(do.call(rbind, dat_st))
  names(dat_st) <- c("Pred_St", "Var_St", "Pred_Xb")
  
  results_hour_st[[i]] <- dat_st
  cat(i,"\n")
}

st_dat <- do.call(rbind,results_hour_st)
st_dat$hour <- rep(12:(12+nh-1), each=nrow(pred_block))

save(st_dat, file=paste0("st_map_preds", dayidx,".Rda"))
