AIC.oglmx<-function(object, ..., k=2){
  # 2*number of estimatated parameters - 2*log likelihood
  value<-k*length(object$coefficients)-2*logLik(object)
  return(value)
}
