## Li, Z., & Ryu, E. (2026). R functions to compute effect size measures in longitudinal growth models. [Manuscript under review]. 
## The r2lmer() function calculates R^2 statistics from models estimated using 
## lmer() in the 'lme4' package (tested on version 2.0-6).
## 
## Requirements:
## - Requires the 'lme4' package.
## - Designed for models with a random intercept (with or without random slopes).
## 
## Arguments:
## - model:   An 'lmer' model object.
## - effectf: Vector of fixed effect terms for which R^2 is computed. 
##            Set to NULL (effectf <- NULL) if no fixed effects are selected.
## - effectr: Vector of random effect terms for which R^2 is computed. 
##            Set to NULL (effectr <- NULL) if no random effects are selected; 
##            use "(Intercept)" to specify the random intercept.


r2lmer<-function(model,effectf,effectr) {
  ## check model fitting issue
  lmer_issues <- function(fit, print = TRUE) {
    stopifnot(inherits(fit, "merMod"))
    oi   <- fit@optinfo
    n    <- length(attr(fit@pp$X, "col.dropped"))
    code <- oi$conv$opt
    
    warningout <- unique(c(
      if (n) sprintf("fixed-effect model matrix is rank deficient so dropping %d %s",
                     n, if (n == 1) "column / coefficient" else "columns / coefficients"),
      tryCatch({ lme4:::checkScaleX(lme4::getME(fit, "X")); NULL },
               warning = conditionMessage),
      unlist(oi$warnings),
      if (!is.null(code) && code != 0)
        paste0("convergence code ", code, " from ", oi$optimizer,
               if (length(oi$message) && nzchar(oi$message)) paste0(": ", oi$message)),
      if (isSingular(fit)) "boundary (singular) fit: see help('isSingular')",
      unlist(oi$conv$lme4$messages)
    ))
    
    if (print && length(warningout)) message(paste(warningout, collapse = "\n"))
    invisible(warningout)
  }
  ## data
  data<-as.matrix(getME(model,"X"))
  datar<-as.matrix(getME(model,"mmList")[[1]])
  ## fixed effects in the model
  xnames<-colnames(model@pp$X)
  ## random effects in the model
  randomeff<-VarCorr(model)
  cid<-names(randomeff)
  xrnames<-colnames(randomeff[[cid]])
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
      sigmaxr<-var(datar[,xrnames])
    }
  } else {
    sigmaxr<-cov(datar[,xrnames])
  }
  ## mean vector of predictors associated with level-2 random effects
  if(length(xrnames)==1) {
    mu<-mean(datar[,xrnames])
  } else {
    mu<-colMeans(datar[,xrnames])
  }
  ## parameter estimates: fixed effect coefficients (gamma)
  gammaall<-as.matrix(getME(model,"beta"))
  gamma<-as.matrix(gammaall[-1])
  ## parameter estimates: covariance matrix of level-2 random effects (tau)
  allrandom<-as.data.frame(VarCorr(model))
  allrandom$var2<-ifelse(is.na(allrandom$var2)==TRUE,allrandom$var1,
                         allrandom$var2)
  tau<-matrix(nrow=length(xrnames),ncol=length(xrnames))
  for (i in 1:nrow(tau)) {
    for (j in 1:ncol(tau)) {
      if (i>=j) {
        tau[i,j]<-allrandom[allrandom$grp==names(model@cnms) & allrandom$var1==xrnames[j] & allrandom$var2==xrnames[i],"vcov"]
      }
      else {
        tau[i,j]<-allrandom[allrandom$grp==names(model@cnms) & allrandom$var1==xrnames[i] & allrandom$var2==xrnames[j],"vcov"]
      }
    }
  }
  rownames(tau)<-xrnames
  colnames(tau)<-xrnames
  ## parameter estimates: variance of level-1 residuals
  esigma2<-allrandom[allrandom$grp=="Residual","vcov"]
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
  names(output)<-c("All Fixed Effects in the model","All Random Effects in the model",
                   "Variance partitioning",
                   "Selected Fixed Effects","Selected Random Effects",
                   "R2 for selected effects",
                   "Warning for correlated terms (Fixed)",
                   "Warning for correlated terms (Random)")
  print(output)
  # lmer warnings and messages
  fitissue<-lmer_issues(model,print=FALSE)
  if (length(fitissue))
    message("\nlmerWarning: lmer() reported the following when fitting this model. ",
            "Check the model before interpreting the results.\n",
            paste(fitissue, collapse = "\n"))
  invisible(output)
}
