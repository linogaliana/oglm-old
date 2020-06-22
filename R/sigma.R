#' Recover residual standard deviation from ordered discrete model
#'
#' Extract the estimated standard deviation of the errors,
#'   the “residual standard deviation” (misnomed also “residual standard error”).
#'
#' This function transforms the linear model for
#'   the standard deviation into the standard deviation
#'
#' @param object A \code{oglmx} model
#' @param ... Additional arguments. Consider in particular adding
#'  `newdata` to parameters
#' @return Residual estimated standard deviation in vector form. With an
#'  homoskedastic model, all values are equal
#' @importFrom stats sigma
#' @export

sigma.oglmx <- function(object, ...){

  args <- list(...)


  if (!inherits(object, "oglmx")) stop("'object' is not a 'oglmx' object")

  if ('newdata' %in% names(args)){
    newdata <- args[['newdata']]
  } else{
    newdata <- NULL
  }


  # TERMS THAT ARE USED FOR VARIANCE COMPUTATION
  # ---------------------------------------------------

  if (is.null(newdata)){
    Z <- object$modelframes$Z
  } else{
    Z <- variance_model(object, newdata = newdata)
  }

  # delta
  # ----------------------------------------------------

  delta <- object$allparams$delta


  # Expression for variance computation
  # ---------------------------------------------------

  if (delta != 0){
    z <- Z %*% delta
  } else{
    z <- 0
  }

  # Return sigma = g(delta*z)
  # --------------------------------

  sigma <- eval(object$sdmodel)

  return(sigma)
}
