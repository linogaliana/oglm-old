probit.reg<-function(formula,data,start=NULL,weights = NULL,beta=NULL,analhessian=TRUE,na.action,savemodelframe=FALSE,robust=FALSE){
  call<-match.call()
  m<-match(c("formula","data","subset","weights", "na.action", "offset"),names(call),0)
  mf<-call[c(1L,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]]<- quote(stats::model.frame)
  mf<-eval(mf,parent.frame())
  weights<-as.vector(model.weights(mf))
  termsMEAN<-terms(mf)
  X<-model.matrix(formula,mf)
  Y<-model.response(mf,"numeric")
  KeepY<-!is.na(Y)
  KeepX<-!apply(X,1,.IsNARow)
  if (!is.null(weights)){
    weights<-weights[KeepX & KeepY]
    robust<-TRUE
    X<-X[!is.na(weights), ,drop=FALSE]
    Y<-Y[!is.na(weights)]
    weights<-weights[!is.na(weights)]
    weighted<-TRUE
    weights<-weights*length(weights)/sum(weights)
  } else {
    weights<-rep_len(1,length(Y))
    weighted<-FALSE
  }
  NoVarModData<-data.frame(Y,weights)

  No.Obs<-length(Y)
  termsMODEL<-list(termsMEAN)
  formulaMODEL<-list(formula)
  # collect variable means and check which variables are binary.
  XVarMeans<-apply(X,2,mean)
  XVarBinary<-apply(X,2,.checkbinary)

  Z<-as.matrix(rep(1,nrow(X)),ncol=1)
  ZVarMeans<-NULL
  ZVarBinary<-NULL
  Heteroskedastic<-FALSE
  BothMeanVar=NULL

  checkoutcomes<-.checkoutcomes(Y,Force=FALSE,binary=TRUE)
  listoutcomes<-checkoutcomes[[1]]
  no.outcomes<-2

  output<-oglmx.fit(Y,X,Z,w=weights,link = "probit",beta=beta,delta=0,threshparam=0,analhessian=analhessian,robustmatrix=robust,start=start)
  output<-append(output,list(call=call,terms=termsMODEL,formula=formulaMODEL,NoVarModData=NoVarModData,Hetero=Heteroskedastic,BothEq=BothMeanVar,varMeans=list(XVarMeans,ZVarMeans),varBinary=list(XVarBinary,ZVarBinary)))
  class(output)<-"oglmx"
  output
}

logit.reg<-function(formula,data,start=NULL,weights=NULL,beta=NULL,analhessian=TRUE,na.action,savemodelframe=FALSE,robust=FALSE){
  call<-match.call()
  m<-match(c("formula","data","subset","weights", "na.action", "offset"),names(call),0)
  mf<-call[c(1L,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]]<- quote(stats::model.frame)
  mf<-eval(mf,parent.frame())
  weights<-as.vector(model.weights(mf))
  termsMEAN<-terms(mf)
  X<-model.matrix(formula,mf)
  Y<-model.response(mf,"numeric")
  KeepY<-!is.na(Y)
  KeepX<-!apply(X,1,.IsNARow)
  if (!is.null(weights)){
    weights<-weights[KeepX & KeepY]
    robust<-TRUE
    X<-X[!is.na(weights), ,drop=FALSE]
    Y<-Y[!is.na(weights)]
    weights<-weights[!is.na(weights)]
    weighted<-TRUE
    weights<-weights*length(weights)/sum(weights)
  } else {
    weights<-rep_len(1,length(Y))
    weighted<-FALSE
  }
  NoVarModData<-data.frame(Y,weights)

  No.Obs<-length(Y)
  termsMODEL<-list(termsMEAN)
  formulaMODEL<-list(formula)
  # collect variable means and check which variables are binary.
  XVarMeans<-apply(X,2,mean)
  XVarBinary<-apply(X,2,.checkbinary)

  Z<-as.matrix(rep(1,nrow(X)),ncol=1)
  ZVarMeans<-NULL
  ZVarBinary<-NULL
  Heteroskedastic<-FALSE
  BothMeanVar=NULL

  checkoutcomes<-.checkoutcomes(Y,Force=FALSE,binary=TRUE)
  listoutcomes<-checkoutcomes[[1]]
  no.outcomes<-2

  output<-oglmx.fit(Y,X,Z,w=weights,link = "logit",beta=beta,delta=0,threshparam=0,analhessian=analhessian,robustmatrix=robust,start=start)
  output<-append(output,list(call=call,terms=termsMODEL,formula=formulaMODEL,NoVarModData=NoVarModData,Hetero=Heteroskedastic,BothEq=BothMeanVar,varMeans=list(XVarMeans,ZVarMeans),varBinary=list(XVarBinary,ZVarBinary)))
  class(output)<-"oglmx"
  output
}

oprobit.reg<-function(formula,data,start=NULL,weights=NULL,beta=NULL,threshparam=NULL,analhessian=TRUE,na.action,savemodelframe=FALSE,robust=FALSE,Force=FALSE){
  call<-match.call()
  m<-match(c("formula","data","subset","weights", "na.action", "offset"),names(call),0)
  mf<-call[c(1L,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]]<- quote(stats::model.frame)
  mf<-eval(mf,parent.frame())
  weights<-as.vector(model.weights(mf))
  termsMEAN<-terms(mf)
  X<-model.matrix(formula,mf)
  Y<-model.response(mf,"numeric")
  KeepY<-!is.na(Y)
  KeepX<-!apply(X,1,.IsNARow)
  if (!is.null(weights)){
    weights<-weights[KeepX & KeepY]
    robust<-TRUE
    X<-X[!is.na(weights), ,drop=FALSE]
    Y<-Y[!is.na(weights)]
    weights<-weights[!is.na(weights)]
    weighted<-TRUE
    weights<-weights*length(weights)/sum(weights)
  } else {
    weights<-rep_len(1,length(Y))
    weighted<-FALSE
  }
  NoVarModData<-data.frame(Y,weights)

  No.Obs<-length(Y)
  # calculate log-likelihood if probabilities are derived from proportions in sample
  # used for calculation of McFadden's R^2
  # BaselineLL<-.BaseLL(Y,weights)
  termsMODEL<-list(termsMEAN)
  formulaMODEL<-list(formula)
  # If specified that there is no constant in the model remove it from the data frame.
  savecolnames<-colnames(X)[colnames(X)!="(Intercept)"]
  X<-X[,colnames(X)!="(Intercept)",drop=FALSE]
  colnames(X)<-savecolnames
  # collect variable means and check which variables are binary.
  XVarMeans<-apply(X,2,mean)
  XVarBinary<-apply(X,2,.checkbinary)

  Z<-as.matrix(rep(1,nrow(X)),ncol=1)
  ZVarMeans<-NULL
  ZVarBinary<-NULL
  Heteroskedastic<-FALSE

  BothMeanVar=NULL
  checkoutcomes<-.checkoutcomes(Y,Force=Force)
  listoutcomes<-checkoutcomes[[1]]
  no.outcomes<-checkoutcomes[[2]]

  output<-oglmx.fit(Y,X,Z,w=weights,link = "probit",beta=beta,delta=0,threshparam=threshparam,analhessian=analhessian,robustmatrix=robust,start=start)
  output<-append(output,list(call=call,terms=termsMODEL,formula=formulaMODEL,NoVarModData=NoVarModData,Hetero=Heteroskedastic,BothEq=BothMeanVar,varMeans=list(XVarMeans,ZVarMeans),varBinary=list(XVarBinary,ZVarBinary)))
  class(output)<-"oglmx"
  output
}

ologit.reg<-function(formula,data,start=NULL,weights=NULL,beta=NULL,threshparam=NULL,analhessian=TRUE,na.action,savemodelframe=FALSE,robust=FALSE,Force=FALSE){
  call<-match.call()
  m<-match(c("formula","data","subset","weights", "na.action", "offset"),names(call),0)
  mf<-call[c(1L,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]]<- quote(stats::model.frame)
  mf<-eval(mf,parent.frame())
  weights<-as.vector(model.weights(mf))
  termsMEAN<-terms(mf)
  X<-model.matrix(formula,mf)
  Y<-model.response(mf,"numeric")
  KeepY<-!is.na(Y)
  KeepX<-!apply(X,1,.IsNARow)
  if (!is.null(weights)){
    weights<-weights[KeepX & KeepY]
    robust<-TRUE
    X<-X[!is.na(weights), ,drop=FALSE]
    Y<-Y[!is.na(weights)]
    weights<-weights[!is.na(weights)]
    weighted<-TRUE
    weights<-weights*length(weights)/sum(weights)
  } else {
    weights<-rep_len(1,length(Y))
    weighted<-FALSE
  }
  NoVarModData<-data.frame(Y,weights)

  No.Obs<-length(Y)
  # calculate log-likelihood if probabilities are derived from proportions in sample
  # used for calculation of McFadden's R^2
  # BaselineLL<-.BaseLL(Y,weights)
  termsMODEL<-list(termsMEAN)
  formulaMODEL<-list(formula)
  # If specified that there is no constant in the model remove it from the data frame.
  savecolnames<-colnames(X)[colnames(X)!="(Intercept)"]
  X<-X[,colnames(X)!="(Intercept)",drop=FALSE]
  colnames(X)<-savecolnames
  # collect variable means and check which variables are binary.
  XVarMeans<-apply(X,2,mean)
  XVarBinary<-apply(X,2,.checkbinary)

  Z<-as.matrix(rep(1,nrow(X)),ncol=1)
  ZVarMeans<-NULL
  ZVarBinary<-NULL
  Heteroskedastic<-FALSE

  BothMeanVar=NULL
  checkoutcomes<-.checkoutcomes(Y,Force=Force)
  listoutcomes<-checkoutcomes[[1]]
  no.outcomes<-checkoutcomes[[2]]

  output<-oglmx.fit(Y,X,Z,w=weights,link = "logit",beta=beta,delta=0,threshparam=threshparam,analhessian=analhessian,robustmatrix=robust,start=start)
  output<-append(output,list(call=call,terms=termsMODEL,formula=formulaMODEL,NoVarModData=NoVarModData,Hetero=Heteroskedastic,BothEq=BothMeanVar,varMeans=list(XVarMeans,ZVarMeans),varBinary=list(XVarBinary,ZVarBinary)))
  class(output)<-"oglmx"
  output
}
