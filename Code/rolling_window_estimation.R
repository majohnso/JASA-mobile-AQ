#####################################################################
#################     Moving Window Estimation     ##################
#####################################################################
for(rep in 1:sum(unique(daysint)>lag)){
  bidx = tblockidx = unique(daysint)[which(unique(daysint)>lag)][rep]
  est_data <- subset(df_block,daysint%in%(tblockidx-1:lag))
  
  byday <- split(est_data, est_data$day)
  only_after_2hrs <- lapply(byday, function(u) subset(u, t > min(t) + 2/24))
  
  est_data2 <- do.call(rbind, only_after_2hrs)
  idx <- which(est_data$idx %in% est_data2$idx)
  
  locstime <- est_data[,c("Longitude","Latitude","t")]
  locstimeX <- est_data[,c("Longitude","Latitude","t",paste0("PC",1:7))]
  linmod = lm(formula(paste("Y_block_med ~ ",paste0(regressorname,collapse ="+"))),data=est_data,na.action = na.exclude)
  est_data$resid = residuals(linmod)
  y <- est_data$resid
  total_variance = var(residuals(linmod))
  
  # space-time model
  startparms <- c(9,0.01,0.01) 
  if(type!=1){
    fit_st_w <- optim(log(startparms), fun_to_max_vecchia.st, var = total_variance,y=y, idx=idx, locstime=locstime, h=h[j], m=m[jj], ncore=ncore, control=list(maxit=1e3, trace=10))
  }else{
    fit_st_w <- optim(log(startparms), fun_to_max_vecchia.st_smple, var = total_variance,y=y, idx=idx, locstime=locstime, h=h[j], m=m[jj], ncore=ncore, nsize=ns, control=list(maxit=1e3, trace=10))
  }
  # space-time-x model
  startparms <- c(9,0.01,0.01,2) 
  if(type!=1){
    fit_stx_w <- optim(log(startparms), fun_to_max_vecchia.stx, var = total_variance,y=y, idx=idx, locstime=locstimeX, h=h[j], m=m[jj], ncore=ncore, control=list(maxit=1e3, trace=10))
  }else{
    fit_stx_w <- optim(log(startparms), fun_to_max_vecchia.stx_smple, var = total_variance,y=y, idx=idx, locstime=locstimeX, h=h[j], m=m[jj], ncore=ncore, nsize=ns, control=list(maxit=1e3, trace=10))
  }
  #### get predictions
  datpred <- subset(df_block,day %in% (wk_block[tblockidx]+0:6))
  
  byday <- split(datpred, floor(datpred$t))
  only_after_1hrs <- lapply(byday, function(u) subset(u, t > min(t) + 1/24))
  
  datpred2 <- do.call(rbind, only_after_1hrs)
  id_pred <- which(datpred$t %in% datpred2$t)
  
  parms_st = exp(fit_st_w$par)
  parms_st[1] = parms_st[1]/(parms_st[1]+1)
  parms_st = c(parms_st[1]*total_variance,parms_st[2],parms_st[3],1,(1-parms_st[1])*total_variance)
  
  parms_stx = exp(fit_stx_w$par)
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
    
    cat("\r",hh)
  }
  
  mspes = rbind(mspe_5min, mspe_15min, mspe_60min)
  cors =  rbind(cors_5min, cors_15min, cors_60min)
  
  save(linmod,pred_5min,pred_15min,pred_60min,mspes,cors,fit_st_w,fit_stx_w,bidx,
       file=paste0("composite_rolling",lag,"wks_lambda_type",type,"h",h[j],"_m",m[jj],"_week",bidx,".Rda"))
}

