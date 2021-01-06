logLik.oglmx<-function(object, ...){
  value<-object$loglikelihood[1]
  attr(value,"df")<-length(object$coefficients)
  return(value)
}

.BaseLL<-function(object){
  data<-object$NoVarModData
  if (ncol(data)==2){
    BaseLL<-as.numeric(logLik(oglmx(Y~1,data=data,weights=data$weights)))
  } else {
    BaseLL<-as.numeric(logLik(oglmx(Y~1,data=data)))
  }
  return(BaseLL)
}


logLik.summary.oglmx<-function(object, ...){
  value<-object$loglikelihood[1]
  attr(value,"df")<-length(object$coefficients)
  return(value)
}
