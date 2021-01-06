oglmx_wrap_null <- function(object, ...){

  if (!is.null(object$call$constantSD)){
    constantSD <- object$call$constantSD
  } else if ('constantSD' %in% names(list(...))){
    constantSD <- list(...)[['constantSD']]
  } else{
    constantSD <- TRUE
  }
  if (!is.null(object$call$SameModelMEANSD)){
    SameModelMEANSD <- object$call$SameModelMEANSD
  } else if ('SameModelMEANSD' %in% names(list(...))){
    SameModelMEANSD <- list(...)[['SameModelMEANSD']]
  } else{
    SameModelMEANSD <- FALSE
  }


  return(
    oglmx(object$NoVarModData$Y~1,
        data=object$NoVarModData,
        weights=object$NoVarModData$weights,
        link = object$link,
        sdmodel = object$sdmodel,
        constantMEAN = "(Intercept)" %in% names(object$coefficients),
        SameModelMEANSD = SameModelMEANSD,
        constantSD = constantSD,
        threshparam = object$allparams$threshparam)
  )

}

