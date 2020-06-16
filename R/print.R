#' @export

print.oglmx <- function(x, digits = max(3L, getOption("digits") - 3L),
                        ...){

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  if (length(coef(x))){

    cat("Coefficients")

    if (is.character(co <- x$contrasts))
      cat("  [contrasts: ", apply(cbind(names(co), co),
                                  1L, paste, collapse = "="), "]")
    cat(":\n")
    print.default(format(x$coefficients, digits = digits),
                  print.gap = 2, quote = FALSE)
  } else{
    cat("No coefficients\n\n")
  }

  cat("\n")

  if (nzchar(mess <- naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")
  cat("Log likelihood by observation:\t   ",
      format(signif(x$loglikelihood/attr(x$loglikelihood, "No.Obs"),
                    digits)),
      "\nAIC:\t", format(signif(AIC(x), digits)))
  cat("\n")
  invisible(x)

}
