logLik.oglmx<-function(object, ...){
  value<-object$loglikelihood[1]
  attr(value,"df")<-length(object$coefficients)
  return(value)
}

.BaseLL<-function(object, ...){
  BaseLL<-as.numeric(logLik(oglmx_wrap_null(object, ...)))
  return(BaseLL)
}


logLik.summary.oglmx<-function(object, ...){
  value<-object$loglikelihood[1]
  attr(value,"df")<-length(object$coefficients)
  return(value)
}


