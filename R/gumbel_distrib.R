#' Gumbel distribution functions
#'
#' Density, distribution function, quantile function
#' and random generation for the gumbel distribtion. Extracted
#' from \link[MASS]{polr} source code
#' @inheritParams stats::rnorm
#' @param loc,scale Location and scale parameters.
#'
#' @return `dgumbel` gives the density, `pgumbel` gives the distribution function,
#' and `rgumbel` generates random deviates.
#'
#' @export

pgumbel <- function(q, loc = 0, scale = 1, lower.tail = TRUE)
{
  q <- (q - loc)/scale
  p <- exp(-exp(-q))
  if (!lower.tail) 1 - p else p
}

#' @rdname pgumbel
#' @export

dgumbel <- function (x, loc = 0, scale = 1, log = FALSE)
{
  x <- (x - loc)/scale
  d <- log(1/scale) - x - exp(-x)
  if (!log) exp(d) else d
}

#' @rdname pgumbel
#' @export

pGumbel <- function(q, loc = 0, scale = 1, lower.tail = TRUE)
{
  q <- (q - loc)/scale
  p <- exp(-exp(q))
  if (lower.tail) 1 - p else p
}

#' @rdname pgumbel
#' @export

dGumbel <- function (x, loc = 0, scale = 1, log = FALSE)
{
  x <- -(x - loc)/scale
  d <- log(1/scale) - x - exp(-x)
  if (!log) exp(d) else d
}

#' @rdname pgumbel
#' @importFrom stats rexp
#' @export

rgumbel <- function(n, loc = 0, scale = 1) loc - scale*log(rexp(n))

#' @rdname pgumbel
#' @export

rGumbel <- function(n, loc = 0, scale = 1) scale*log(rexp(n)) - loc
