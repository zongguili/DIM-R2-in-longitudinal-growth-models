## Li, Z., & Ryu, E. (2026). R functions to compute effect size measures in longitudinal growth models. [Manuscript under review].
## The r2lmeRmat() function calculates R^2 statistics from models estimated using 
## lme() in the 'nlme' package (tested on version 3.1-171).
## The function supports the following specifications for 
## the level-1 residual covariance structure (R matrix): Constant, Diag, CS, CSH, AR1, ARH1. 
## - Constant: constant variance
## - Diag: diagonal structure with heterogeneous variances
## - CS: Compound symmetry
## - CSH: Compound symmetry with heterogeneous variances 
## - AR1: First-order autocorrelation
## - ARH1: First-order autocorrelation with heterogeneous variances
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
## r2lmeRmat() provides an additional check for singular fits (threshold: 1e-4).

r2lmeRmat<-function(model,effectf,effectr) {
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
  ## level-1 residual covariance structure
  corSatt<-attributes(model$modelStruct$corStruct)$class[1]
  varSatt<-attributes(model$modelStruct$varStruct)$class[1]
  if (is.null(corSatt)==TRUE && is.null(varSatt)==TRUE) {
    Rmat<-"Constant"
    } else if (is.null(corSatt)==TRUE && varSatt=="varIdent") {
        Rmat<-"Diag"
    } else if (corSatt=="corCompSymm" && is.null(varSatt)==TRUE) {
        Rmat<-"CS"
    } else if (corSatt=="corCompSymm" && varSatt=="varIdent") {
        Rmat<-"CSH"
    } else if (corSatt=="corAR1" && is.null(varSatt)==TRUE) {
        Rmat<-"AR1"
    } else if (corSatt=="corAR1" && varSatt=="varIdent") {
        Rmat<-"ARH1"
    } else {
        Rmat<-NULL
    }
  ## parameter estimates: covariance matrix of level-1 residuals
  if (Rmat=="Constant"){
    esigma2<-(model$sigma)^2
  } else if (Rmat=="Diag") {
    formula<-as.character(attributes(model$modelStruct$varStruct)$formula)[2]
    timevar<-sub(".*\\| ","",formula)
    sigmasqtt<-(coef(model$modelStruct$varStruct,unconstrained=F,allCoef=T)*(model$sigma))^2
    esigma<-diag(sigmasqtt)
  } else if (Rmat=="CS") {
    formula<-as.character(attributes(model$modelStruct$corStruct)$formula)[2]
    timevar<-sub(" \\|.*","",formula)
    rho<-coef(model$modelStruct$corStruct,unconstrained=F)
    csmat<-matrix(rho,nrow=nlevels(factor(data[[timevar]])),ncol=nlevels(factor(data[[timevar]])))
    diag(csmat)<-1
    esigma<-csmat*(model$sigma)^2
  } else if (Rmat=="CSH") {
    formula<-as.character(attributes(model$modelStruct$varStruct)$formula)[2]
    timevar<-sub(".*\\| ","",formula)
    sigmatt<-coef(model$modelStruct$varStruct,unconstrained=F,allCoef=T)*(model$sigma)
    rho<-coef(model$modelStruct$corStruct,unconstrained=F)
    csmat<-matrix(rho,nrow=nlevels(factor(data[[timevar]])),ncol=nlevels(factor(data[[timevar]])))
    diag(csmat)<-1
    esigma<-diag(sigmatt)%*%csmat%*%diag(sigmatt)
  } else if (Rmat=="AR1") {
    formula<-as.character(attributes(model$modelStruct$corStruct)$formula)[2]
    timevar<-sub(" \\|.*","",formula)
    phi<-coef(model$modelStruct$corStruct,unconstrained=F)
    ar1cor<-function(n,phi) {
      exponent<-abs(matrix(1:n-1,nrow=n,ncol=n,byrow=TRUE)-(1:n-1))
      phi^exponent
    }
    armat<-ar1cor(nlevels(factor(data[[timevar]])),phi)
    esigma<-armat*(model$sigma)^2
  } else if (Rmat=="ARH1") {
    formula<-as.character(attributes(model$modelStruct$varStruct)$formula)[2]
    timevar<-sub(".*\\| ","",formula)
    sigmatt<-coef(model$modelStruct$varStruct,unconstrained=F,allCoef=T)*(model$sigma)
    phi<-coef(model$modelStruct$corStruct,unconstrained=F)
    ar1cor<-function(n,phi) {
      exponent<-abs(matrix(1:n-1,nrow=n,ncol=n,byrow=TRUE)-(1:n-1))
      phi^exponent
    }
    armat<-ar1cor(nlevels(factor(data[[timevar]])),phi)
    esigma<-diag(sigmatt)%*%armat%*%diag(sigmatt)
  } else {
    esigma<-NULL
  }
  ## level-1 residual variance unweighted and weighted
  if (Rmat=="Constant"){
    residvar<-esigma2
    residvarw<-esigma2
  } else if (is.null(Rmat)==FALSE) {
    nt<-table(data[,timevar])
    nt<-as.vector(nt)
    residvar<-sum(diag(esigma))*(1/nlevels(factor(data[[timevar]])))
    ntw<-nt/sum(nt)
    esigmaw<-diag(sqrt(ntw))%*%esigma%*%diag(sqrt(ntw))
    residvarw<-sum(diag(esigmaw))
  }
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
  vartotal<-t(gamma)%*%sigmax%*%gamma+sum(diag(tau%*%sigmaxr))+t(mu)%*%tau%*%mu+residvar
  vartotalw<-t(gamma)%*%sigmax%*%gamma+sum(diag(tau%*%sigmaxr))+t(mu)%*%tau%*%mu+residvarw
  varf<-t(gamma)%*%sigmax%*%gamma
  varr<-sum(diag(tau%*%sigmaxr))+t(mu)%*%tau%*%mu
  vareff<-(t(DIMf%*%gamma)%*%sigmax%*%(DIMf%*%gamma)+(sum(diag((DIMr%*%tau%*%DIMr)%*%sigmaxr))+t(mu)%*%(DIMr%*%tau%*%DIMr)%*%mu))
  ## R^2 
  R2f<-varf/vartotal
  R2r<-varr/vartotal
  R2e<-residvar/vartotal
  R2<-vareff/vartotal
  R2fw<-varf/vartotalw
  R2rw<-varr/vartotalw
  R2ew<-residvarw/vartotalw
  R2w<-vareff/vartotalw
  value1<-c(varf,varr,residvar,vartotal,R2f,R2r,R2e)
  names(value1)<-c("AllFixed","AllRandom","Residual","TotalVar","R2AllFixed","R2AllRandom","R2Residual")
  value2<-c(vartotal,vareff,R2)
  names(value2)<-c("TotalVar","EffectVar","R2")
  value1w<-c(varf,varr,residvarw,vartotalw,R2fw,R2rw,R2ew)
  names(value1w)<-c("AllFixed","AllRandom","Residual","TotalVar","R2AllFixed","R2AllRandom","R2Residual")
  value2w<-c(vartotalw,vareff,R2w)
  names(value2w)<-c("TotalVar","EffectVar","R2")
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
  if (Rmat=="Constant"|Rmat=="CS"|Rmat=="AR1") {
    outputuw<-list(xnames,xrnames,Rmat,value1,effectf,effectr,value2,warnf,warnr)
    names(outputuw)<-c("All Fixed Effects in the model",
                       "All Random Effects in the model",
                       "Level-1 residual covariance structure",
                       "Variance partitioning",
                       "Selected Fixed Effects","Selected Random Effects",
                       "R2 for selected effects",
                       "Warning for correlated terms (Fixed)",
                       "Warning for correlated terms (Random)")
    output<-outputuw
  } else if (Rmat=="Diag"|Rmat=="CSH"|Rmat=="ARH1") {
    outputw<-list(xnames,xrnames,Rmat,nt,value1,value1w,effectf,effectr,value2,value2w,warnf,warnr)
    names(outputw)<-c("All Fixed Effects in the model",
                      "All Random Effects in the model",
                      "Level-1 residual covariance structure",
                      "Number of observations at each time point",
                      "Variance partitioning",
                      "Variance partitioning (weighted residual variance)",
                      "Selected Fixed Effects","Selected Random Effects",
                      "R2 for selected effects",
                      "R2 for selected effects (weighted residual variance)",
                      "Warning for correlated terms (Fixed)",
                      "Warning for correlated terms (Random)")
    output<-outputw
    }
  print(output)
  # lme warnings and messages; singular fit check
  fitissue<-lme_issues(model,print=FALSE)
  if (length(fitissue))
    message("lmeWarning: lme() reported the following when fitting this model. Check the model before interpreting the results.\n",
            paste(fitissue, collapse = "\n"))
  invisible(output)
}
