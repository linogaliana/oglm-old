vcov.oglmx<-function(object,tol=1e-20,...){
  if (is.null(object$BHHHhessian)){
    vcov<-qr.solve(-object$hessian,tol=tol)
  } else {
    vcov<- qr.solve(object$hessian,tol=tol)%*%(object$BHHHhessian*(attr(object$loglikelihood,"No.Obs")/(attr(object$loglikelihood,"No.Obs")-1)))%*%qr.solve(object$hessian,tol=tol)
  }
  colnames(vcov)<-rownames(vcov)<-names(object$coefficients)
  return(vcov)
}
