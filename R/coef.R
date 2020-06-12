coef.oglmx<-function(object, ...){
  coefnames<-names(object$coefficients)
  output<-as.vector(object$coefficients)
  names(output)<-coefnames
  return(output)
}

coef.summary.oglmx<-function(object, ...){
  coefnames<-names(object$coefficients)
  output<-as.vector(object$coefficients)
  names(output)<-coefnames
  return(output)
}
