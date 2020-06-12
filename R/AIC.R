#' Akaike's An Information Criterion for `oglm` objects
#'
#'
#' @inheritParams stats::AIC
#'
#' @description Compute `AIC` for ordered and interval
#'  regression objects. See [stats::AIC()] for the rest
#'
#' @details See [stats::AIC()]
#' @seealso [stats::AIC()],
#'  [stats::logLik()], [stats::nobs()]

AIC.oglmx<-function(object, ..., k=2){
  # 2*number of estimatated parameters - 2*log likelihood
  value<-k*length(object$coefficients)-2*logLik(object)
  return(value)
}
