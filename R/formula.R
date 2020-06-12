#' Obtain model formula for an `oglm` object
#'
#' Given an object of class `oglm` the function describes
#' the estimated model via an expression of class [stats::formula()].
#' The function serves to provide a name of a model to the [lmtest::lrtest()]
#'
#' @param x object of class `oglmx`
#' @param ... Additional arguments, currently ignored.
#' @return an object of class `formula`.


formula.oglmx<-function(x, ...){
  # extract the formula for an oglmx object
  # for use to apply a model name in lrtest
  if (is.null(x$formula[[2]])){
    value<-x$formula[[1]]
  } else {
    # collect the names from the terms output
    # from the mean equation, include response term
    meannames<-names(attr(terms(x)[[1]],"dataClasses"))
    varnames<-attr(terms(x)[[2]],"term.labels")
    textoutput<-paste(meannames[1],"~",meannames[2])
    if (length(meannames)>2){
      for (j in 3:length(meannames)){
        textoutput<-paste(textoutput,"+",meannames[j])
      }
    }
    textoutput<-paste(textoutput,"|",varnames[1])
    if (length(varnames)>1){
      for (j in 2:length(varnames)){
        textoutput<-paste(textoutput,"+",varnames[j])
      }
    }
    value<-formula(textoutput)
  }
  return(value)
}
