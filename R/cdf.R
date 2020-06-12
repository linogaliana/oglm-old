.cdf.func<-function(link){
  if (link=="probit"){
    value <- function(p){pnorm(p)}
  } else if (link=="logit"){
    value <- function(p){plogis(p)}
  } else if (link=="cauchit"){
    value <- function(p){pcauchy(p)}
  } else if (link=="loglog"){
    value <-function(p){exp(-exp(-p))}
  } else if (link=="cloglog"){
    value <- function(p){1-exp(-exp(p))}
  } else {
    stop("Specified link not available.")
  }
  value
}

.pdf.func<-function(link){
  if (link=="probit"){
    value <- function(p){dnorm(p)}
  } else if (link=="logit"){
    value <- function(p){dlogis(p)}
  } else if (link=="cauchit"){
    value <- function(p){dcauchy(p)}
  } else if (link=="loglog"){
    value <-function(p){exp(-exp(-p))*exp(-p)}
  } else if (link=="cloglog"){
    value <- function(p){exp(-exp(p))*exp(p)}
  } else {
    stop("Specified link not available.")
  }
  value
}

.Dpdf.func<-function(link){
  if (link=="probit"){
    value <- function(p){-p*dnorm(p)}
  } else if (link=="logit"){
    value <- function(p){dlogis(p)*(1-2*plogis(p))}
  } else if (link=="cauchit"){
    value <- function(p){-2*p*dcauchy(p)/(1+p^2)}
  } else if (link=="loglog"){
    value <-function(p){exp(-exp(-p))*(exp(-p)*(exp(-p)-1))}
  } else if (link=="cloglog"){
    value <- function(p){exp(-exp(p))*exp(p)*(1-exp(p))}
  } else {
    stop("Specified link not available.")
  }
  value
}

.DDpdf.func<-function(link){
  if (link=="probit"){
    value <- function(p){(p^2-1)*dnorm(p)}
  } else if (link=="logit"){
    value <- function(p){dlogis(p)*(1-2*plogis(p))^2-2*dlogis(p)^2}
  } else if (link=="cauchit"){
    value <- function(p){(dcauchy(p)/(p^2+1)^4)*(6*p^4+8*p^2-2)}
  } else if (link=="loglog"){
    value <- function(p){exp(-exp(-p))*(exp(-3*p)-3*exp(-2*p)+exp(-p))}
  } else if (link=="cloglog"){
    value <- function(p){exp(-exp(p))*(exp(3*p)-3*exp(2*p)+exp(p))}
  } else {
    stop("Specified link not available.")
  }
  value
}
