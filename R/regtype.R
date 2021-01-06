.regtype.oglmx<-function(object){

  if (ncol(object$NoVarModData)==2){

    if (sum(object$NoVarModData$weights==1)!=nrow(object$NoVarModData)){ # this is needed for backwards compatibility, for saved data before v3
      NotWeight<-FALSE
    } else {
      NotWeight<-TRUE
    }
  } else {
    NotWeight<-TRUE

  }

  Zero<-ifelse(NotWeight,"","Weighted ")
  First<-ifelse(object$Hetero,"Heteroskedastic ","")
  Second<-ifelse(object$NOutcomes>2,"Ordered ","")

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
