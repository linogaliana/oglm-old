.regtype.oglmx<-function(object){
  if (sum(object$NoVarModData$weights==1)!=nrow(object$NoVarModData)){
    Zero<-"Weighted "
  } else {
    Zero<-""
  }
  if (object$Hetero){
    First<-"Heteroskedastic "
  } else {
    First<-""
  }
  if (object$NOutcomes>2){
    Second<-"Ordered "
  } else {
    Second<-""
  }
  if (object$link=="logit"){
    Third<-"Logit "
  } else if (object$link=="probit"){
    Third<-"Probit "
  } else if (object$link=="cloglog"){
    Third<-"CLogLog "
  } else if (object$link=="loglog"){
    Third<-"LogLog "
  } else if (object$link=="cauchit"){
    Third<-"Cauchit "
  }
  Fourth<-"Regression"
  value<-paste(Zero,First,Second,Third,Fourth,sep="")
  return(value)
}
