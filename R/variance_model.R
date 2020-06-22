variance_model <- function(object, newdata = NULL,
                           ...){

  if (!inherits(object, "oglmx")) stop("'object' should be a 'oglmx' object")

  # WHEN newdata is NULL, we use the initial object
  if (is.null(newdata)) newdata <- object$modelframes$Z


  if (!is.null(object$formula$sdeq)) {
    # if (object$constantSD) {
    #   Z <- newdata
    #   Zint <- match("(Intercept)", colnames(Z), nomatch = 0L)
    #   if (Zint > 0L) {
    #     Z <- Z[, -Zint, drop = FALSE]
    #   }
    # } else{
    Z <- newdata
    # }
    Z <- model.matrix(object$formula$sdeq, data.frame(Z))
  }
  else {
    Z <- matrix(rep(1, nrow(newdata)), ncol = 1)
  }


  return(Z)
}
