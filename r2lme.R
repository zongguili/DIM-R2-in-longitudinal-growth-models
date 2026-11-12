## Li, Z., & Ryu, E. (2026). R functions to compute effect size measures in longitudinal growth models. [Manuscript under review].
## The r2lme() function calculates R^2 statistics from models estimated using 
## lme() in the 'nlme' package (tested on version 3.1-171).
## 
## Requirements:
## - Requires the `nlme` package.
## - Designed for models with a random intercept (with or without random slopes).
## 
## Arguments:
## - model:   An 'lme' model object.
## - effectf: Vector of fixed effect terms for which R^2 is computed. 
##            Set to NULL (effectf <- NULL) if no fixed effects are selected.
## - effectr: Vector of random effect terms for which R^2 is computed. 
##            Set to NULL (effectr <- NULL) if no random effects are selected; 
##            use "(Intercept)" to specify the random intercept.
##
## r2lme() provides an additional check for singular fits (threshold: 1e-4).


r2lme<-function(model,effectf,effectr) {
  ## check model fitting issue
  lme_issues <- function(model, print = TRUE, tol = 1e-4) {
    rel_sd <- sapply(as.matrix(model$modelStruct$reStruct), function(m)
      sqrt(max(0, min(eigen(m, symmetric = TRUE, only.values = TRUE)$values))))
    
    warningout <- c(
      if (is.character(model$apVar)) paste(model$apVar),
      if (any(rel_sd < tol))
        paste("Additional check (not reported by lme()): possible singular fit;",
              "a random-effect variance is near zero or a correlation is near +/-1 (see VarCorr())")
    )
    
    if (print && length(warningout)) message(paste(warningout, collapse = "\n"))
    invisible(warningout)
  }
  ## data
  if (is.null(model$na.action)==TRUE) {
    data<-model$data
  } else {
    data<-model$data
    data<-data[!(rownames(data)) %in% names(model$na.action),]
  }
  data[,"(Intercept)"]<-1
  xnames<-names(model$coefficients$fixed)
  for (xname in xnames) {
    if (grepl(":",xname)==TRUE){
      x1<-sub("\\:.*", "",xname)
      x2<-sub(".*:", "",xname)
      intname<-paste(x1,":",x2,sep="")
      data[,intname]<-data[,x1]*data[,x2]
    } else {
    }
  }
  cid<-names(model$coefficients$random)
  xrnames<-colnames(model$coefficients$random[[cid]])
  ## covariance matrix of predictors (sigmax, sigmaxr) 
  ## sigmax
  if(length(xnames)==1) {
    if(xnames=="(Intercept)") {
      sigmax<-0
    } else {
      sigmax<-var(data[,xnames])
    }
  } else if (length(xnames)==2 & "(Intercept)" %in% xnames){
    sigmax<-as.matrix(var(data[,xnames[xnames!="(Intercept)"]]))
    rownames(sigmax)<-xnames[xnames!="(Intercept)"]
    colnames(sigmax)<-xnames[xnames!="(Intercept)"]
  } else {
    sigmax<-cov(data[,xnames[xnames!="(Intercept)"]])
  }
  ## sigmaxr
  if(length(xrnames)==1) {
    if(xrnames=="(Intercept)") {
      sigmaxr<-0
    } else {
      sigmaxr<-var(data[,xrnames])
    }
  } else{
    sigmaxr<-cov(data[,xrnames])
  }
  ## mean vector of predictors associated with level-2 random effects
  if(length(xrnames)==1){
    mu<-mean(data[,xrnames])
  } else {
    mu<-colMeans(data[,xrnames])
  }
  ## parameter estimates: fixed effect coefficients (gamma)
  gammaall<-model$coefficients$fixed
  gamma<-as.matrix(gammaall[names(gammaall)!="(Intercept)"])
  ## parameter estimates: covariance matrix of level-2 random effects (tau)
  tau<-getVarCov(model)
  ## parameter estimates: variance of level-1 residuals
  esigma2<-(model$sigma)^2
  ## Dummy indicator matrix for fixed effects 
  if(length(gamma)==0) {
    gamma<-0
    DIMf<-0
  } else {
    DIMf<-matrix(0,nrow=length(gamma),ncol=1)
    rownames(DIMf)<-xnames[xnames!="(Intercept)"]
    DIMf[is.element(rownames(DIMf),effectf)]<-1
    if (nrow(DIMf)>1) {
      DIMf<-diag(as.vector(DIMf))
    } else {
      rownames(DIMf)<-NULL
    }
  }
  ## Dummy indicator matrix for random effects 
  DIMr<-matrix(0,nrow=length(xrnames),ncol=1)
  rownames(DIMr)<-xrnames
  DIMr[is.element(rownames(DIMr),effectr)]<-1
  if (nrow(DIMr)>1) {
    DIMr<-diag(as.vector(DIMr))
  } else {
    rownames(DIMr)<-NULL
  }
  ## variance 
  vartotal<-t(gamma)%*%sigmax%*%gamma+sum(diag(tau%*%sigmaxr))+t(mu)%*%tau%*%mu+esigma2
  varf<-t(gamma)%*%sigmax%*%gamma
  varr<-sum(diag(tau%*%sigmaxr))+t(mu)%*%tau%*%mu
  vareff<-(t(DIMf%*%gamma)%*%sigmax%*%(DIMf%*%gamma)+(sum(diag((DIMr%*%tau%*%DIMr)%*%sigmaxr))+t(mu)%*%(DIMr%*%tau%*%DIMr)%*%mu))
  ## R^2 
  R2f<-varf/vartotal
  R2r<-varr/vartotal
  R2e<-esigma2/vartotal
  R2<-vareff/vartotal
  value1<-c(varf,varr,esigma2,vartotal,R2f,R2r,R2e)
  names(value1)<-c("AllFixed","AllRandom","Residual","TotalVar","R2AllFixed","R2AllRandom","R2Residual")
  value2<-c(vartotal,vareff,R2)
  names(value2)<-c("TotalVar","EffectVar","R2")
  ## warning for correlated terms
  fxrnames<-colnames(sigmax)[(colnames(sigmax) %in% effectf)]
  fxcnames<-colnames(sigmax)[!(colnames(sigmax) %in% effectf)]
  fx<-sigmax[fxrnames,fxcnames]
  rx<-tau[rownames(tau)==effectr,colnames(tau)!=effectr]
  if (!is.null(fx) && !all(fx==0)) {
    warnf<-"The selected fixed effect term(s) are correlated with other terms."
  } else {
    warnf<-NA 
  }
  if (!is.null(rx) && !all(rx==0)) {
    warnr<-"The selected random effect term(s) are correlated with other terms."
  } else {
    warnr<-NA 
  }
  output<-list(xnames,xrnames,value1,effectf,effectr,value2,warnf,warnr)
  names(output)<-c("All Fixed Effects in the model",
                   "All Random Effects in the model",
                   "Variance partitioning",
                   "Selected Fixed Effects",
                   "Selected Random Effects",
                   "R2 for selected effects",
                   "Warning for correlated terms (Fixed)",
                   "Warning for correlated terms (Random)")
  print(output)
  # lme warnings and messages; singular fit check
  fitissue<-lme_issues(model,print=FALSE)
  if (length(fitissue))
    message("lmeWarning: lme() reported the following when fitting this model. Check the model before interpreting the results.\n",
            paste(fitissue, collapse = "\n"))
  invisible(output)
}
