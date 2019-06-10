##### spatial estimation and kriging functions

s_covfn.gp_cpp <- function(parms,st1,st2){
  if(identical(st1,st2)){
    if(nrow(st1)==1){
      s1 <- matrix(unlist(st1[,1:2]),nr=1)
      s2 <- matrix(unlist(st2[,1:2]),nr=1)
    }else{
      s1 <- as.matrix(st1[,1:2])
      s2 <- as.matrix(st2[,1:2])
    }
    s <- matern_s_cpp(parms[1], parms[2], 0.5, s1, s2, sq=TRUE)
    diag(s) <- parms[1] + parms[3]
  }else{
    s1 <- as.matrix(st1[,1:2])
    s2 <- as.matrix(st2[,1:2])
    if(nrow(st1)==1){
      s1 <- matrix(unlist(st1[,1:2]),nr=1)
    }
    if(nrow(st2)==1){
      s2 <- matrix(unlist(st2[,1:2]),nr=1)
    }
    s <- matern_s_cpp(parms[1], parms[2], 0.5, s1, s2, sq=FALSE)
  }
  return(s)
}

logl.gp.mh.s <- function(i, y, locstime, parms, m){ # m are minutes 
  # h is lag to first neighbor, m is interval past lag h
  ub <- locstime[i,3] - 1/(60*60*24)
  lb <- locstime[i,3] - m/60/24
  int <- dplyr::between(locstime[,3], lb, ub)
  
  # int <- (i-1):(i-m)
  
  if(sum(int) > 0){
    condat <- as.matrix(cbind(locstime[int,],y[int]))
    if(sum(int)==1){
      condat <- matrix(condat, nr=1)
    }
    preddat <- matrix(c(locstime[i,],y[i]),nr=1)
    
    # s1 pred loc, s2 cond loc
    s12 <- s_covfn.gp_cpp(parms,preddat,condat)
    s2 <- s_covfn.gp_cpp(parms,condat,condat)
    s1 <- s_covfn.gp_cpp(parms,preddat,preddat) - s12%*%solve(s2,t(s12))
    y2 <- unlist(condat[,ncol(condat)])
    y1 <- unlist(preddat[,ncol(preddat)])
    yhat <- y1 - s12%*%solve(s2,y2)
    # logl <- -n/2*log(2*pi)-1/2*determinant(s1,log=T)$mod-1/2*t(yhat)%*%solve(s1,yhat)
    logl <- -1/2*log(2*pi)-1/2*log(s1)-1/2*yhat^2/s1
    return(logl)
  }else{
    return(NA)
  }
}

krig.mh.s <- function(i, nbs, datpred, tr_dat, cov_names, parms_s, linmod){
  # h is lag to first neighbor, m is interval past lag h
  # ub <- datpred[i,"t"] - h/60/24
  # lb <- datpred[i,"t"] - (h+m)/60/24
  
  int <- nbs[i,]
  
  # Y2|Y1
  X1 = tr_dat[int,cov_names] 
  X2 = data.frame(datpred[i,cov_names])
  Y1 = tr_dat[int,"Y_block_med"]
  # Y2 = datpred[i,"Y_block_med"]
  s1 = tr_dat[int,c("Longitude","Latitude")]
  s2 = data.frame(datpred[i,c("Longitude","Latitude")])
  n = nrow(s1)
  
  if(sum(int)==1){
    X1 <- data.frame(X1)
    s1 <- data.frame(s1)
  }
  
  Y1hat = predict(linmod, newdata = cbind(s1,X1))
  Y2hat = predict(linmod, newdata = cbind(s2,X2))
  
  pr_s    <- krig_STX_cpp(parms_s,as.matrix(s2[,1:2]), rep(1,nrow(s2)), matrix(1, nr=nrow(s2), nc=1),as.matrix(s1[,1:2]),rep(1,nrow(s1)),matrix(1, nr=nrow(s1), nc=1),Y1,Y2hat, Y1hat)
  as.numeric(c(pr_s$pred, pr_s$var, Y2hat))
}

vecchia_logl.gp.mh.s<-function(y, idx, parms, locstime, m, ncore){
  # N <- length(y)
  # idx <- seq(m+h+1,N,by=1)
  logl_all <- mclapply(idx, function(x) logl.gp.mh.s(x, y, locstime, parms, m),mc.cores=ncore)
  # logl_all <- foreach(j=1:length(idx),.noexport=c("matern_st_cpp")) %dopar% {logl.gp.mh(idx[j], y, locstime,parms,h,m)}
  sum(unlist(logl_all),na.rm=TRUE)
}

fun_to_max_vecchia.s <- function(logparms, var, y, idx, locstime, m, ncore){
  parms <- rep(0,2)
  parms[2] <- exp(logparms[2])
  parms[1] = exp(logparms[1])/(exp(logparms[1])+1)
  parms = c(parms[1]*var,parms[2],(1-parms[1])*var)
  -vecchia_logl.gp.mh.s(y, idx, parms, locstime, m, ncore) + 0.01*sum(parms^2)
}

logl.gp.mh.s_smple <- function(i, y, locstime, parms, m, nsize=800){ # h and m are minutes 
  # h is lag to first neighbor, m is interval past lag h
  ub <- locstime[i,3]
  lb <- locstime[i,3] - m/60/24
  int <- dplyr::between(locstime[,3], lb, ub)
  if(sum(int)!=0){
    int <- sample(which(int), min(nsize, sum(int)))
  }
  int <- (1:nrow(locstime)) %in% int
  
  if(sum(int) > 0){
    condat <- as.matrix(cbind(locstime[int,],y[int]))
    if(sum(int)==1){
      condat <- matrix(condat, nr=1)
    }
    preddat <- matrix(c(locstime[i,],y[i]),nr=1)
    
    # s1 pred loc, s2 cond loc
    s12 <- s_covfn.gp_cpp(parms,preddat,condat)
    s2 <- s_covfn.gp_cpp(parms,condat,condat)
    s1 <- s_covfn.gp_cpp(parms,preddat,preddat) - s12%*%solve(s2,t(s12))
    y2 <- unlist(condat[,ncol(condat)])
    y1 <- unlist(preddat[,ncol(preddat)])
    yhat <- y1 - s12%*%solve(s2,y2)
    # logl <- -n/2*log(2*pi)-1/2*determinant(s1,log=T)$mod-1/2*t(yhat)%*%solve(s1,yhat)
    logl <- -1/2*log(2*pi)-1/2*log(s1)-1/2*yhat^2/s1
    return(logl)
  }else{
    return(NA)
  }
}

vecchia_logl.gp.mh.s_smple<-function(y, idx, parms, locstime, h, m, ncore, nsize=800){
  # N <- length(y)
  # idx <- seq(m+h+1,N,by=1)
  logl_all <- mclapply(idx, function(x) logl.gp.mh.s_smple(x, y, locstime, parms, m, nsize=nsize),mc.cores=ncore)
  # logl_all <- foreach(j=1:length(idx),.noexport=c("matern_st_cpp")) %dopar% {logl.gp.mh(idx[j], y, locstime,parms,h,m)}
  sum(unlist(logl_all),na.rm=TRUE)
}

fun_to_max_vecchia.s_smple <- function(logparms, var, y, idx, locstime, h, m, ncore, nsize=800){
  parms <- rep(0,2)
  parms[2] <- exp(logparms[2])
  parms[1] = exp(logparms[1])/(exp(logparms[1])+1)
  parms = c(parms[1]*var,parms[2],(1-parms[1])*var)
  -vecchia_logl.gp.mh.s_smple(y, idx, parms, locstime, h, m, ncore, nsize=nsize)
}

##### space-time estimation and kriging functions
st_covfn.gp_cpp <- function(parms,st1,st2){
  if(identical(st1,st2)){
    if(nrow(st1)==1){
      s1 <- matrix(unlist(st1[,1:2]),nr=1)
      s2 <- matrix(unlist(st2[,1:2]),nr=1)
    }else{
      s1 <- as.matrix(st1[,1:2])
      s2 <- as.matrix(st2[,1:2])
    }
    s <- matern_st_cpp(parms[1], parms[2], parms[3], 0.5, as.numeric(st1[,3]), as.numeric(st2[,3]), s1, s2, sq=TRUE)
    diag(s) <- parms[1] + parms[4]
  }else{
    s1 <- as.matrix(st1[,1:2])
    s2 <- as.matrix(st2[,1:2])
    if(nrow(st1)==1){
      s1 <- matrix(unlist(st1[,1:2]),nr=1)
    }
    if(nrow(st2)==1){
      s2 <- matrix(unlist(st2[,1:2]),nr=1)
    }
    s <- matern_st_cpp(parms[1], parms[2], parms[3], 0.5, as.numeric(st1[,3]), as.numeric(st2[,3]), s1, s2, sq=FALSE)
  }
  return(s)
}

logl.gp.mh.st <- function(i, y, locstime, parms, h, m){ # h and m are minutes 
  # h is lag to first neighbor, m is interval past lag h
  ub <- locstime[i,3] - h/60/24
  lb <- locstime[i,3] - (h+m)/60/24
  int <- dplyr::between(locstime[,3], lb, ub)
  
  if(sum(int) > 0){
    condat <- as.matrix(cbind(locstime[int,],y[int]))
    if(sum(int)==1){
      condat <- matrix(condat, nr=1)
    }
    preddat <- matrix(c(locstime[i,],y[i]),nr=1)
    
    # s1 pred loc, s2 cond loc
    s12 <- st_covfn.gp_cpp(parms,preddat,condat)
    s2 <- st_covfn.gp_cpp(parms,condat,condat)
    s1 <- st_covfn.gp_cpp(parms,preddat,preddat) - s12%*%solve(s2,t(s12))
    y2 <- unlist(condat[,ncol(condat)])
    y1 <- unlist(preddat[,ncol(preddat)])
    yhat <- y1 - s12%*%solve(s2,y2)
    # logl <- -n/2*log(2*pi)-1/2*determinant(s1,log=T)$mod-1/2*t(yhat)%*%solve(s1,yhat)
    logl <- -1/2*log(2*pi)-1/2*log(s1)-1/2*yhat^2/s1
    return(logl)
  }else{
    return(NA)
  }
}

krig.mh.st <- function(i, datpred, parms, h, m, cov_names, linmod){
  # h is lag to first neighbor, m is interval past lag h
  ub <- datpred[i,"t"] - h/60/24
  lb <- datpred[i,"t"] - (h+m)/60/24
  
  int <- dplyr::between(datpred[,"t"], lb, ub)
  
  # Y2|Y1
  X1 = datpred[int,cov_names] 
  X2 = data.frame(datpred[i,cov_names])
  Y1 = datpred[int,"Y_block_med"]
  Y2 = datpred[i,"Y_block_med"]
  s1 = datpred[int,c("Longitude","Latitude","t")]
  s2 = data.frame(datpred[i,c("Longitude","Latitude","t")])
  n = nrow(s1)
  
  if(sum(int)==1){
    X1 <- data.frame(X1)
    s1 <- data.frame(s1)
  }
  
  Y1hat = predict(linmod, newdata = cbind(s1,X1))
  Y2hat = predict(linmod, newdata = cbind(s2,X2))
  
  # pred, cond
  pr_st   <- krig_STX_cpp(parms,as.matrix(s2[,1:2]), s2$t, matrix(1, nr=nrow(s2), nc=1),as.matrix(s1[,1:2]),s1$t,matrix(1, nr=nrow(s1), nc=1),Y1,Y2hat, Y1hat)
  as.numeric(c(pr_st$pred, pr_st$var, Y2hat))
}

vecchia_logl.gp.mh.st<-function(y, idx, parms, locstime, h, m, ncore){
  # N <- length(y)
  # idx <- seq(m+h+1,N,by=1)
  logl_all <- mclapply(idx, function(x) logl.gp.mh.st(x, y, locstime, parms, h, m),mc.cores=ncore)
  # logl_all <- foreach(j=1:length(idx),.noexport=c("matern_st_cpp")) %dopar% {logl.gp.mh(idx[j], y, locstime,parms,h,m)}
  sum(unlist(logl_all),na.rm=TRUE)
}

fun_to_max_vecchia.st <- function(logparms, var, y, idx, locstime, h, m, ncore){
  parms <- rep(0,3)
  parms[2:3] <- exp(logparms[2:3])
  parms[1] = exp(logparms[1])/(exp(logparms[1])+1)
  parms = c(parms[1]*var,parms[2],parms[3],(1-parms[1])*var)
  -vecchia_logl.gp.mh.st(y, idx, parms, locstime, h, m, ncore) + 0.01*sum(parms^2)
}

logl.gp.mh.st_smple <- function(i, y, locstime, parms, h, m, nsize=800){ # h and m are minutes 
  # h is lag to first neighbor, m is interval past lag h
  ub <- locstime[i,3] - h/60/24
  lb <- locstime[i,3] - (h+m)/60/24
  int <- dplyr::between(locstime[,3], lb, ub)
  if(sum(int)!=0){
    int <- sample(which(int), min(nsize, sum(int)))
  }
  int <- (1:nrow(locstime)) %in% int
  
  if(sum(int) > 0){
    condat <- as.matrix(cbind(locstime[int,],y[int]))
    if(sum(int)==1){
      condat <- matrix(condat, nr=1)
    }
    preddat <- matrix(c(locstime[i,],y[i]),nr=1)
    
    # s1 pred loc, s2 cond loc
    s12 <- st_covfn.gp_cpp(parms,preddat,condat)
    s2 <- st_covfn.gp_cpp(parms,condat,condat)
    s1 <- st_covfn.gp_cpp(parms,preddat,preddat) - s12%*%solve(s2,t(s12))
    y2 <- unlist(condat[,ncol(condat)])
    y1 <- unlist(preddat[,ncol(preddat)])
    yhat <- y1 - s12%*%solve(s2,y2)
    # logl <- -n/2*log(2*pi)-1/2*determinant(s1,log=T)$mod-1/2*t(yhat)%*%solve(s1,yhat)
    logl <- -1/2*log(2*pi)-1/2*log(s1)-1/2*yhat^2/s1
    return(logl)
  }else{
    return(NA)
  }
}

krig.mh.st_smple <- function(i, datpred, parms, h, m, cov_names, linmod, nsize=800){
  # h is lag to first neighbor, m is interval past lag h
  ub <- datpred[i,"t"] - h/60/24
  lb <- datpred[i,"t"] - (h+m)/60/24
  
  int <- dplyr::between(datpred[,"t"], lb, ub)
  if(sum(int)!=0){
    int <- sample(which(int), min(nsize, sum(int)))
  }
  int <- (1:nrow(datpred)) %in% int
  
  # Y2|Y1
  X1 = datpred[int,cov_names] 
  X2 = data.frame(datpred[i,cov_names])
  Y1 = datpred[int,"Y_block_med"]
  Y2 = datpred[i,"Y_block_med"]
  s1 = datpred[int,c("Longitude","Latitude","t")]
  s2 = data.frame(datpred[i,c("Longitude","Latitude","t")])
  n = nrow(s1)
  
  if(sum(int)==1){
    X1 <- data.frame(X1)
    s1 <- data.frame(s1)
  }
  
  Y1hat = predict(linmod, newdata = cbind(s1,X1))
  Y2hat = predict(linmod, newdata = cbind(s2,X2))
  
  # pred, cond
  pr_st   <- krig_STX_cpp(parms,as.matrix(s2[,1:2]), s2$t, matrix(1, nr=nrow(s2), nc=1),as.matrix(s1[,1:2]),s1$t,matrix(1, nr=nrow(s1), nc=1),Y1,Y2hat, Y1hat)
  as.numeric(c(pr_st$pred, pr_st$var, Y2hat))
}

vecchia_logl.gp.mh.st_smple<-function(y, idx, parms, locstime, h, m, ncore, nsize=800){
  # N <- length(y)
  # idx <- seq(m+h+1,N,by=1)
  logl_all <- mclapply(idx, function(x) logl.gp.mh.st_smple(x, y, locstime, parms, h, m, nsize=nsize),mc.cores=ncore)
  # logl_all <- foreach(j=1:length(idx),.noexport=c("matern_st_cpp")) %dopar% {logl.gp.mh(idx[j], y, locstime,parms,h,m)}
  sum(unlist(logl_all),na.rm=TRUE)
}

fun_to_max_vecchia.st_smple <- function(logparms, var, y, idx, locstime, h, m, ncore, nsize=800){
  parms <- rep(0,3)
  parms[2:3] <- exp(logparms[2:3])
  parms[1] = exp(logparms[1])/(exp(logparms[1])+1)
  parms = c(parms[1]*var,parms[2],parms[3],(1-parms[1])*var)
  -vecchia_logl.gp.mh.st_smple(y, idx, parms, locstime, h, m, ncore, nsize=nsize)
}

krig.mh.st_carA <- function(i, datpred, parms, cov_names, linmod){
  # h is lag to first neighbor, m is interval past lag h
  datpredB <- datpred[datpred$Car==2,]
  datpredA <- datpred[datpred$Car==1,]
  
  intB <- dplyr::between(datpredB[,"t"], i - 1/24, i + 1/24 - 1e-7)
  intA <- dplyr::between(datpredA[,"t"], i, i + 1/24 - 1e-7)
  
  if(sum(intA) > 0){
    # Y2|Y1
    X1 = datpredB[intB,cov_names] 
    X2 = data.frame(datpredA[intA,cov_names])
    Y1 = datpredB[intB,"Y_block_med"]
    Y2 = datpredA[intA,"Y_block_med"]
    s1 = datpredB[intB,c("Longitude","Latitude","t")]
    s2 = data.frame(datpredA[intA,c("Longitude","Latitude","t")])
    n = nrow(s1)
    
    if(sum(intB)==1){
      X1 <- data.frame(X1)
      s1 <- data.frame(s1)
    }
    
    Y1hat = predict(linmod, newdata = cbind(s1,X1))
    Y2hat = predict(linmod, newdata = cbind(s2,X2))
    
    # pred, cond
    pr_st   <- krig_STX_cpp(parms,as.matrix(s2[,1:2]), s2$t, matrix(1, nr=nrow(s2), nc=1),as.matrix(s1[,1:2]),s1$t,matrix(1, nr=nrow(s1), nc=1),Y1,Y2hat, Y1hat)
    return(data.frame(pred_st = pr_st$pred, var_st = pr_st$var, Y2hat))
  }else{
    return(NULL)
  }
}

krig.mh.st_carA_smple <- function(i, datpred, parms, cov_names, linmod, nsize=800){
  # h is lag to first neighbor, m is interval past lag h
  datpredB <- datpred[datpred$Car==2,]
  datpredA <- datpred[datpred$Car==1,]
  
  intB <- dplyr::between(datpredB[,"t"], i - 1/24, i + 1/24 - 1e-7)
  if(sum(intB)!=0){
    intB <- sample(which(intB), min(nsize, sum(intB)))
  }
  intB <- (1:nrow(datpredB)) %in% intB
  
  intA <- dplyr::between(datpredA[,"t"], i, i + 1/24 - 1e-7)
  
  if(sum(intA) > 0){
    # Y2|Y1
    X1 = datpredB[intB,cov_names] 
    X2 = data.frame(datpredA[intA,cov_names])
    Y1 = datpredB[intB,"Y_block_med"]
    Y2 = datpredA[intA,"Y_block_med"]
    s1 = datpredB[intB,c("Longitude","Latitude","t")]
    s2 = data.frame(datpredA[intA,c("Longitude","Latitude","t")])
    n = nrow(s1)
    
    if(sum(intB)==1){
      X1 <- data.frame(X1)
      s1 <- data.frame(s1)
    }
    
    Y1hat = predict(linmod, newdata = cbind(s1,X1))
    Y2hat = predict(linmod, newdata = cbind(s2,X2))
    
    # pred, cond
    pr_st   <- krig_STX_cpp(parms,as.matrix(s2[,1:2]), s2$t, matrix(1, nr=nrow(s2), nc=1),as.matrix(s1[,1:2]),s1$t,matrix(1, nr=nrow(s1), nc=1),Y1,Y2hat, Y1hat)
    return(data.frame(pred_st = pr_st$pred, var_st = pr_st$var, Y2hat))
  }else{
    return(NULL)
  }
}


##### space-time-x estimation and kriging functions
stx_covfn.gp_cpp <- function(parms,st1,st2){
  if(identical(st1,st2)){
    if(nrow(st1)==1){
      s1 <- matrix(unlist(st1[,1:2]),nr=1)
      s2 <- matrix(unlist(st2[,1:2]),nr=1)
      
      x1 <- matrix(unlist(st1[,4:ncol(st1)]),nr=1)
      x2 <- matrix(unlist(st2[,4:ncol(st1)]),nr=1)
    }else{
      s1 <- as.matrix(st1[,1:2])
      s2 <- as.matrix(st2[,1:2])
      
      x1 <- as.matrix(st1[,4:ncol(st1)])
      x2 <- as.matrix(st2[,4:ncol(st1)])
    }
    s <- matern_stx_cpp(parms[1], parms[2], parms[3], parms[4], 0.5, as.numeric(st1[,3]), as.numeric(st2[,3]), x1, x2, s1, s2, sq=TRUE)
    diag(s) <- parms[1] + parms[5]
  }else{
    s1 <- as.matrix(st1[,1:2])
    s2 <- as.matrix(st2[,1:2])
    
    x1 <- as.matrix(st1[,4:ncol(st1)])
    x2 <- as.matrix(st2[,4:ncol(st1)])
    if(nrow(st1)==1){
      
      s1 <- matrix(unlist(st1[,1:2]),nr=1)
      x1 <- matrix(unlist(st1[,4:ncol(st1)]),nr=1)
    }
    if(nrow(st2)==1){
      s2 <- matrix(unlist(st2[,1:2]),nr=1)
      x2 <- matrix(unlist(st2[,4:ncol(st1)]),nr=1)
    }
    s <- matern_stx_cpp(parms[1], parms[2], parms[3], parms[4], 0.5, as.numeric(st1[,3]), as.numeric(st2[,3]), x1, x2, s1, s2, sq=FALSE)
  }
  return(s)
}

logl.gp.mh.stx <- function(i, y, locstime, parms, h, m){ # h and m are minutes 
  # h is lag to first neighbor, m is interval past lag h
  ub <- locstime[i,3] - h/60/24
  lb <- locstime[i,3] - (h+m)/60/24
  int <- dplyr::between(locstime[,3], lb, ub)
  
  if(sum(int) > 0){
    condat <- as.matrix(cbind(locstime[int,],y[int]))
    if(sum(int)==1){
      condat <- matrix(condat, nr=1)
    }
    preddat <- matrix(c(locstime[i,],y[i]),nr=1)
    
    # s1 pred loc, s2 cond loc
    s12 <- stx_covfn.gp_cpp(parms,preddat,condat)
    s2 <- stx_covfn.gp_cpp(parms,condat,condat)
    s1 <- stx_covfn.gp_cpp(parms,preddat,preddat) - s12%*%solve(s2,t(s12))
    y2 <- unlist(condat[,ncol(condat)])
    y1 <- unlist(preddat[,ncol(preddat)])
    yhat <- y1 - s12%*%solve(s2,y2)
    # logl <- -n/2*log(2*pi)-1/2*determinant(s1,log=T)$mod-1/2*t(yhat)%*%solve(s1,yhat)
    logl <- -1/2*log(2*pi)-1/2*log(s1)-1/2*yhat^2/s1
    return(logl)
  }else{
    return(NA)
  }
}

krig.mh.stx <- function(i, datpred, parms, h, m, cov_names, linmod){
  # h is lag to first neighbor, m is interval past lag h
  ub <- datpred[i,"t"] - h/60/24
  lb <- datpred[i,"t"] - (h+m)/60/24
  
  int <- dplyr::between(datpred[,"t"], lb, ub)
  
  # Y2|Y1
  X1 = datpred[int,cov_names] 
  X2 = data.frame(datpred[i,cov_names])
  Y1 = datpred[int,"Y_block_med"]
  Y2 = datpred[i,"Y_block_med"]
  s1 = datpred[int,c("Longitude","Latitude","t")]
  s2 = data.frame(datpred[i,c("Longitude","Latitude","t")])
  n = nrow(s1)
  
  if(sum(int)==1){
    X1 <- data.frame(X1)
    s1 <- data.frame(s1)
  }
  
  Y1hat = predict(linmod, newdata = cbind(s1,X1))
  Y2hat = predict(linmod, newdata = cbind(s2,X2))
  
  pr_stx   <- krig_STX_cpp(parms,as.matrix(s2[,1:2]), s2$t, as.matrix(X2[,1:7]),as.matrix(s1[,1:2]),s1$t,as.matrix(X1[,1:7]),Y1,Y2hat, Y1hat)
  as.numeric(c(pr_stx$pred, pr_stx$var, Y2hat))
}

vecchia_logl.gp.mh.stx<-function(y, idx, parms, locstime, h, m, ncore){
  # N <- length(y)
  # idx <- seq(m+h+1,N,by=1)
  logl_all <- mclapply(idx, function(x) logl.gp.mh.stx(x, y, locstime, parms, h, m),mc.cores=ncore)
  # logl_all <- foreach(j=1:length(idx),.noexport=c("matern_st_cpp")) %dopar% {logl.gp.mh(idx[j], y, locstime,parms,h,m)}
  sum(unlist(logl_all),na.rm=TRUE)
}

fun_to_max_vecchia.stx <- function(logparms, var, y, idx, locstime, h, m, ncore){
  parms <- rep(0,4)
  parms[2:4] <- exp(logparms[2:4])
  parms[1] = exp(logparms[1])/(exp(logparms[1])+1)
  parms = c(parms[1]*var,parms[2:4],(1-parms[1])*var)
  -vecchia_logl.gp.mh.stx(y, idx, parms, locstime, h, m, ncore) + 0.01*sum(parms^2)
}

logl.gp.mh.stx_smple <- function(i, y, locstime, parms, h, m, nsize=800){ # h and m are minutes 
  # h is lag to first neighbor, m is interval past lag h
  ub <- locstime[i,3] - h/60/24
  lb <- locstime[i,3] - (h+m)/60/24
  int <- dplyr::between(locstime[,3], lb, ub)
  if(sum(int)!=0){
    int <- sample(which(int), min(nsize, sum(int)))
  }
  int <- (1:nrow(locstime)) %in% int
  
  if(sum(int) > 0){
    condat <- as.matrix(cbind(locstime[int,],y[int]))
    if(sum(int)==1){
      condat <- matrix(condat, nr=1)
    }
    preddat <- matrix(c(locstime[i,],y[i]),nr=1)
    
    # s1 pred loc, s2 cond loc
    s12 <- stx_covfn.gp_cpp(parms,preddat,condat)
    s2 <- stx_covfn.gp_cpp(parms,condat,condat)
    s1 <- stx_covfn.gp_cpp(parms,preddat,preddat) - s12%*%solve(s2,t(s12))
    y2 <- unlist(condat[,ncol(condat)])
    y1 <- unlist(preddat[,ncol(preddat)])
    yhat <- y1 - s12%*%solve(s2,y2)
    # logl <- -n/2*log(2*pi)-1/2*determinant(s1,log=T)$mod-1/2*t(yhat)%*%solve(s1,yhat)
    logl <- -1/2*log(2*pi)-1/2*log(s1)-1/2*yhat^2/s1
    return(logl)
  }else{
    return(NA)
  }
}

krig.mh.stx_smple <- function(i, datpred, parms, h, m, cov_names, linmod, nsize=800){
  # h is lag to first neighbor, m is interval past lag h
  ub <- datpred[i,"t"] - h/60/24
  lb <- datpred[i,"t"] - (h+m)/60/24
  
  int <- dplyr::between(datpred[,"t"], lb, ub)
  if(sum(int)!=0){
    int <- sample(which(int), min(nsize, sum(int)))
  }
  int <- (1:nrow(datpred)) %in% int
  
  # Y2|Y1
  X1 = datpred[int,cov_names] 
  X2 = data.frame(datpred[i,cov_names])
  Y1 = datpred[int,"Y_block_med"]
  Y2 = datpred[i,"Y_block_med"]
  s1 = datpred[int,c("Longitude","Latitude","t")]
  s2 = data.frame(datpred[i,c("Longitude","Latitude","t")])
  n = nrow(s1)
  
  if(sum(int)==1){
    X1 <- data.frame(X1)
    s1 <- data.frame(s1)
  }
  
  Y1hat = predict(linmod, newdata = cbind(s1,X1))
  Y2hat = predict(linmod, newdata = cbind(s2,X2))
  
  # pred, cond
  pr_stx   <- krig_STX_cpp(parms,as.matrix(s2[,1:2]), s2$t, as.matrix(X2[,1:7]),as.matrix(s1[,1:2]),s1$t,as.matrix(X1[,1:7]),Y1,Y2hat, Y1hat)
  as.numeric(c(pr_stx$pred, pr_stx$var, Y2hat))
}

vecchia_logl.gp.mh.stx_smple<-function(y, idx, parms, locstime, h, m, ncore, nsize=800){
  # N <- length(y)
  # idx <- seq(m+h+1,N,by=1)
  logl_all <- mclapply(idx, function(x) logl.gp.mh.stx_smple(x, y, locstime, parms, h, m, nsize=nsize),mc.cores=ncore)
  # logl_all <- foreach(j=1:length(idx),.noexport=c("matern_st_cpp")) %dopar% {logl.gp.mh(idx[j], y, locstime,parms,h,m)}
  sum(unlist(logl_all),na.rm=TRUE)
}

fun_to_max_vecchia.stx_smple <- function(logparms, var, y, idx, locstime, h, m, ncore, nsize=800){
  parms <- rep(0,4)
  parms[2:4] <- exp(logparms[2:4])
  parms[1] = exp(logparms[1])/(exp(logparms[1])+1)
  parms = c(parms[1]*var,parms[2:4],(1-parms[1])*var)
  -vecchia_logl.gp.mh.stx_smple(y, idx, parms, locstime, h, m, ncore, nsize=nsize)
}

krig.mh.stx_carA <- function(i, datpred, parms, cov_names, linmod){
  # h is lag to first neighbor, m is interval past lag h
  datpredB <- datpred[datpred$Car==2,]
  datpredA <- datpred[datpred$Car==1,]
  
  intB <- dplyr::between(datpredB[,"t"], i - 1/24, i + 1/24 - 1e-7)
  intA <- dplyr::between(datpredA[,"t"], i, i + 1/24 - 1e-7)
  
  if(sum(intA) > 0){
    # Y2|Y1
    X1 = datpredB[intB,cov_names] 
    X2 = data.frame(datpredA[intA,cov_names])
    Y1 = datpredB[intB,"Y_block_med"]
    Y2 = datpredA[intA,"Y_block_med"]
    s1 = datpredB[intB,c("Longitude","Latitude","t")]
    s2 = data.frame(datpredA[intA,c("Longitude","Latitude","t")])
    n = nrow(s1)
    
    if(sum(intB)==1){
      X1 <- data.frame(X1)
      s1 <- data.frame(s1)
    }
    
    Y1hat = predict(linmod, newdata = cbind(s1,X1))
    Y2hat = predict(linmod, newdata = cbind(s2,X2))
    
    # pred, cond
    pr_stx   <- krig_STX_cpp(parms,as.matrix(s2[,1:2]), s2$t, as.matrix(X2[,1:7]), as.matrix(s1[,1:2]),s1$t, as.matrix(X1[,1:7]),Y1,Y2hat, Y1hat)
    return(data.frame(pred_stx = pr_stx$pred, var_st = pr_stx$var, Y2hat))
  }else{
    return(NULL)
  }
}

krig.mh.stx_carA_smple <- function(i, datpred, parms, cov_names, linmod, nsize=800){
  # h is lag to first neighbor, m is interval past lag h
  datpredB <- datpred[datpred$Car==2,]
  datpredA <- datpred[datpred$Car==1,]
  
  intB <- dplyr::between(datpredB[,"t"], i - 1/24, i + 1/24 - 1e-7)
  if(sum(intB)!=0){
    intB <- sample(which(intB), min(nsize, sum(intB)))
  }
  intB <- (1:nrow(datpredB)) %in% intB
  
  intA <- dplyr::between(datpredA[,"t"], i, i + 1/24 - 1e-7)
  
  if(sum(intA) > 0){
    # Y2|Y1
    X1 = datpredB[intB,cov_names] 
    X2 = data.frame(datpredA[intA,cov_names])
    Y1 = datpredB[intB,"Y_block_med"]
    Y2 = datpredA[intA,"Y_block_med"]
    s1 = datpredB[intB,c("Longitude","Latitude","t")]
    s2 = data.frame(datpredA[intA,c("Longitude","Latitude","t")])
    n = nrow(s1)
    
    if(sum(intB)==1){
      X1 <- data.frame(X1)
      s1 <- data.frame(s1)
    }
    
    Y1hat = predict(linmod, newdata = cbind(s1,X1))
    Y2hat = predict(linmod, newdata = cbind(s2,X2))
    
    # pred, cond
    pr_stx   <- krig_STX_cpp(parms,as.matrix(s2[,1:2]), s2$t, as.matrix(X2[,1:7]), as.matrix(s1[,1:2]),s1$t, as.matrix(X1[,1:7]),Y1,Y2hat, Y1hat)
    return(data.frame(pred_stx = pr_stx$pred, var_st = pr_stx$var, Y2hat))
  }else{
    return(NULL)
  }
}

st_covfn <- function(parms,st1,st2){
  s1 = matrix(NA,nrow=nrow(st1),ncol(st1))
  s2 = matrix(NA,nrow=nrow(st2),ncol(st2))
  
  s1[,1:2]=st1[,1:2]/parms[2]
  s2[,1:2]=st2[,1:2]/parms[2]
  s1[,3]=st1[,3]/parms[3]; s2[,3]=st2[,3]/parms[3]
  d = rdist(s1,s2)
  s <- parms[1]*exp(-d)
  if(identical(st1,st2))  diag(s) = parms[1] + parms[4]
  return(s)
}

mspe <- function(X1,X2,st1,st2,parms){
  s21 <- st_covfn(parms,st2,st1)
  s11 <- st_covfn(parms,st1,st1)
  s22 <- st_covfn(parms,st2,st2)
  var = mean(parms[1]+parms[4] - diag(s21%*%solve(s11,t(s21))))
  return(var)
}

mse <- function(X1,X2,st1,st2,parms){
  s21 <- st_covfn(parms,st2,st1)
  s11 <- st_covfn(parms,st1,st1)
  s22 <- st_covfn(parms,st2,st2)
  var = mean(parms[1]+parms[4] - diag(s21%*%solve(s11,t(s21))))
  return(var)
}

MSE <- function(x) mean(x^2)
