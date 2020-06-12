logLik.oglmx<-function(object, ...){
  value<-object$loglikelihood[1]
  attr(value,"df")<-length(object$coefficients)
  return(value)
}

.BaseLL<-function(object){
  BaseLL<-as.numeric(logLik(oglmx(object$NoVarModData$Y~1,data=object$NoVarModData,weights=object$NoVarModData$weights)))
  return(BaseLL)
}


logLik.summary.oglmx<-function(object, ...){
  value<-object$loglikelihood[1]
  attr(value,"df")<-length(object$coefficients)
  return(value)
}
