McFaddensR2.oglmx<-function(object){
  value<-1-logLik(object)/.BaseLL(object)
  return(value)
}
