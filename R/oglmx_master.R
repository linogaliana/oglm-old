#' Fit Ordered Generalized Linear Model.
#'
#' @description oglmx is used to estimate models for which
#' the outcome variable is discrete and the mean and/or
#' variance of the underlying latent variable can be
#' modelled as a linear combination of explanatory variables.
#' Standard models such as probit, logit, ordered probit
#' and ordered logit are included in the diverse
#' set of models estimated by the function.
#'
#' @param formulaMEAN an object of class \link[stats]{formula}:
#'  a symbolic description of the model used to explain the
#'  mean of the latent variable. The response variable
#'  should be a numeric vector or factor variable
#'  such that the numerical assignments for the
#'  levels of the factor have ordinal meaning.
#' @param formulaSD either `NULL` (homoskedastic model) or an object
#'   of class \link[stats]{formula}: a symbolic description of
#'   the model used to explain the variance of the latent variable.
#' @param data a data frame containing the variables in the model.
#' @param start either `NULL` or a numeric vector specifying
#'   start values for each of the estimated parameters,
#'   passed to the maximisation routine.
#' @param link specifies a link function for the model to be
#'   estimated, accepted values are *"probit"*,
#'   *"logit"*, *"cauchit"*, *"loglog"* and *"cloglog"*
#' @param constantMEAN logical. Should an intercept be
#'   included in the model of the mean of the latent variable?
#'   Can be overwritten and set to `FALSE` using the
#'   `formulaMEAN` argument by writing `0 +`
#'   as the first element of the equation.
#' @param constantSD logical. Should an intercept be included
#'   in the model of the variance of the latent variable?
#'   Can be overwritten and set to `FALSE` using the `formulaSD`
#'   argument by writing `0 +` as the first element of the equation.
#' @param beta `NULL` or numeric vector. Used to prespecify
#'   elements of the parameter vector for the equation of the mean
#'   of the latent variable. Vector should be of length one or of
#'   length equal to the number of explanatory
#'   variables in the mean equation.
#'   If of length one the value is presumed to correspond to the
#'   constant if a constant is included or the first element of the
#'   parameter vector. If of length greater than one then `NA`
#'   should be entered for elements of the vector to be estimated.
#' @param delta `NULL` or numeric vector. Used to prespecify
#'   elements of the parameter vector for the equation of
#'   the variance of the latent variable. Vector should be of
#'   length one or of length equal to the number of explanatory
#'   variables in the variance equation.
#'   If of length one the value is presumed to
#'   correspond to the constant if a constant is included
#'   or the first element of the parameter vector. If
#'   of length greater than one then `NA` should be entered
#'   for elements of the vector to be estimated.
#' @param threshparam `NULL` or numeric vector. Used to prespecify
#'   the threshold parameters of the model. Vector should be
#'   of length equal to the number of outcomes minus one. `NA` should
#'   be entered for threshold parameters to be estimated by the model.
#' @param analhessian logical. Indicates whether the analytic
#'   Hessian should be calculated and used, default is
#'   `TRUE`, if set to `FALSE` a finite-difference
#'   approximation of the Hessian is used.
#' @param sdmodel object of mode *“expression”*. The expression
#'   defines the function that transforms the linear model
#'   for the standard deviation into the standard deviation.
#'   The expression should be written as a function of variable `z`.
#'   The default value is `expression(exp(z))`.
#' @param SameModelMEANSD logical. Indicates whether the matrix used to
#'   model the mean of the latent variable is identical to that used
#'   to model the variance. If `formulaSD=NULL` and `SameModelMEANSD=TRUE` a
#'   model with heteroskedasticity is estimated. If
#' `SameModelMEANSD=FALSE` and `formulaSD==formulaMEAN`
#'   value is overridden. Used to reduce memory requirements when models are identical.
#' @param savemodelframe logical. Indicates whether the
#'   model frame(s) should be saved for future use.
#'   Default is `FALSE`. Should be set to `TRUE` if
#'   intending to estimate Average Marginal Effects.
#' @param Force logical. If set to `FALSE` (the default),
#'   the function stops if the response
#'   variable has more than twenty categories.
#'   Should be changed to `TRUE` if a model with
#'   more than twenty categories is desired.
#' @param robust logical. If set to `TRUE` the outer product or
#'   BHHH estimate of the meat in the sandwich of the
#'   variance-covariance matrix is calculated. If calculated
#'   standard errors will be calculated using the sandwich estimator
#'   by default when calling `summary`.
#' @param optmeth specifies a method for the maximisation of the likelihood
#'   passed to [maxLik::maxLik()]. Default to *NR* (Newton-Raphson)
#' @inheritParams stats::glm
#'
#' @return An object of class "\code{oglmx}" with the following components:
#' \describe{
#'     \item{loglikelihood}{log-likelihood for the estimated model.
#'       Includes as attributes the log-likelihood for the constant
#'       only model and the number of observations.}
#'     \item{link}{link function used in the estimated model.}
#'     \item{no.iterations}{number of iterations of maximisation algorithm.}
#'     \item{coefficients}{named vector of estimated parameters.}
#'     \item{returnCode}{code returned by
#'       the `maxLik` optimisation routine}
#'     \item{call}{the call used to generate the results.}
#'     \item{gradient}{numeric vector, the value of the
#'       gradient of the log-likelihood function at the obtained
#'       parameter vector. Should be approximately equal to zero.}
#'     \item{terms}{two element list. Each element is an object of
#'       type `terms` related to the mean and standard
#'       deviation equation respectively.
#'     }
#'     \item{formula}{two element list. Each element is an
#'       object of type [stats::formula()] related to
#'       the mean and standard deviation
#'       equation respectively.}
#'     \item{NoVarModData}{dataframe. Contains data required to
#'       estimate the no information model used in
#'       calculation of McFadden's R-squared measure.}
#'     \item{hessian}{hessian matrix of the log-likelihood
#'       function evaluated at the obtained parameter vector.}
#'     \item{BHHHhessian}{Either `NULL` if no weights were
#'       included and `robust = FALSE`, or the BHHH estimate.}
#'     \item{Hetero}{logical. If `TRUE` indicates that the
#'       estimated model includes a model for the variance
#'       of the error term, i.e. heteroskedasticity.}
#'     \item{NOutcomes}{the number of distinct outcomes
#'       in the response variable.}
#'     \item{Outcomes}{numeric vector of length equal to `NOutcomes`.
#'       Lists the values of the different outcomes.}
#'     \item{BothEq}{data.frame with either two or three
#'       columns. Lists the names of variables that are in both
#'       the mean and variance equations and their locations
#'       within their respective model frames. Information is required in
#'       the call of `margins.oglmx` to obtain correct marginal effects.}
#'     \item{allparams}{a list containing three numeric vectors,
#'       the vectors contain the parameters from the mean equation,
#'       the variance equation and the threshold parameters respectively.
#'       Includes the prespecified and estimated parameters together.}
#'     \item{varMeans}{a list containing two numeric vectors. The vectors
#'       list the mean values of the variables in the mean and variance
#'       equation respectively. Stored for use in a call of
#'       \code{margins.oglmx} to obtain marginal effects at means.}
#'     \item{varBinary}{a list containing two numeric vectors. The vectors
#'       indicate whether the variables in the mean and variance
#'       equations are binary indicators. Stored for use in a call
#'       of \code{margins.oglmx} to obtain marginal effects at means.}
#'     \item{Est.Parameters}{list containing three logical vectors.
#'       Indicates which parameters in the parameter vectors were estimated.}
#'     \item{modelframes}{If \code{savemodelframe} set to \code{FALSE}
#'       then returns \code{NULL}, otherwise returns a list with two
#'       elements, the model frames for the mean and variance equations.}
#' }
#'
#' @examples \dontrun{
#' # create random sample, three variables, two binary.
#' set.seed(242)
#' n<-250
#' x1<-sample(c(0,1),n,replace=TRUE,prob=c(0.75,0.25))
#' x2<-vector("numeric",n)
#' x2[x1==0]<-sample(c(0,1),n-sum(x1==1),replace=TRUE,prob=c(2/3,1/3))
#' z<-rnorm(n,0.5)
#' # create latent outcome variable
#' latenty<-0.5+1.5*x1-0.5*x2+0.5*z+rnorm(n,sd=exp(0.5*x1-0.5*x2))
#' # observed y has four possible values: -1,0,1,2
#' # threshold values are: -0.5, 0.5, 1.5.
#' y<-vector("numeric",n)
#' y[latenty< -0.5]<- -1
#' y[latenty>= -0.5 & latenty<0.5]<- 0
#' y[latenty>= 0.5 & latenty<1.5]<- 1
#' y[latenty>= 1.5]<- 2
#' dataset<-data.frame(y,x1,x2)
#' # estimate standard ordered probit
#' results.oprob<-oglmx(y ~ x1 + x2 + z, data=dataset,link="probit",constantMEAN=FALSE,
#'                      constantSD=FALSE,delta=0,threshparam=NULL)
#' coef(results.oprob) # extract estimated coefficients
#' summary(results.oprob)
#' # calculate marginal effects at means
#' margins.oglmx(results.oprob)
#' # estimate ordered probit with heteroskedasticity
#' results.oprobhet<-oglmx(y ~ x1 + x2 + z, ~ x1 + x2, data=dataset, link="probit",
#'                         constantMEAN=FALSE, constantSD=FALSE,threshparam=NULL)
#' summary(results.oprobhet)
#' library("lmtest")
#' # likelihood ratio test to compare model with and without heteroskedasticity.
#' lrtest(results.oprob,results.oprobhet)
#' # calculate marginal effects at means.
#' margins.oglmx(results.oprobhet)
#' # scale of parameter values is meaningless. Suppose instead two of the
#' # three threshold values were known, then can include constants in the
#' # mean and standard deviation equation and the scale is meaningful.
#' results.oprobhet1<-oglmx(y ~ x1 + x2 + z, ~ x1 + x2, data=dataset, link="probit",
#'                          constantMEAN=TRUE, constantSD=TRUE,threshparam=c(-0.5,0.5,NA))
#' summary(results.oprobhet1)
#' margins.oglmx(results.oprobhet1)
#' # marginal effects are identical to results.oprobithet, but using the true thresholds
#' # means the estimated parameters are on the same scale as underlying data.
#' # can choose any two of the threshold values and get broadly the same result.
#' results.oprobhet2<-oglmx(y ~ x1 + x2 + z, ~ x1 + x2, data=dataset, link="probit",
#'                          constantMEAN=TRUE, constantSD=TRUE,threshparam=c(-0.5,NA,1.5))
#' summary(results.oprobhet2)
#' margins.oglmx(results.oprobhet2)
#' # marginal effects are again identical. Parameter estimates do change.
#' }
#'
#' @import maxLik
#' @import stats
#' @export

oglmx <- function(formulaMEAN,
                  formulaSD=NULL,
                  data,
                  start=NULL,weights=NULL,
                  link="probit",
                  constantMEAN=TRUE, constantSD=TRUE,
                  beta=NULL,delta=NULL,
                  threshparam=NULL,
                  analhessian=TRUE,
                  sdmodel=expression(exp(z)),
                  optmeth = c("NR", "BFGS", "BFGSR", "BHHH", "SANN", "CG", "NM"),
                  SameModelMEANSD=FALSE,
                  na.action,
                  savemodelframe=TRUE,
                  Force=FALSE,
                  robust=FALSE){

  optmeth <- match.arg(optmeth)
  call <- match.call()

  names(call)[match("formulaMEAN",names(call),0)]<-"formula"
  m<-match(c("formula","data","subset","weights", "na.action", "offset"),names(call),0)
  mf<-call[c(1L,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]]<- quote(stats::model.frame)
  mf<-eval(mf,parent.frame())
  weights<-as.vector(model.weights(mf))
  termsMEAN<-terms(mf)
  X <- model.matrix(as.formula(formulaMEAN),
                  mf)
  Y <- model.response(mf,"numeric")
  Keep<-!is.na(Y) & !apply(X,1,.IsNARow)
  if (!is.null(formulaSD) & !SameModelMEANSD){
    dataframeSD<-model.frame(formulaSD,data=data,na.action=na.pass)
    Z<-model.matrix(as.formula(formulaSD),dataframeSD)
    Keep<-Keep & !apply(Z,1,.IsNARow)
    Z<-Z[Keep, ,drop=FALSE]
    termsSD<-terms(dataframeSD)
    # If specified that there is no constant in the model remove it from the data frame.
    if (!constantSD & ncol(Z)>1){
      savecolnamesZ<-colnames(Z)[colnames(Z)!="(Intercept)"]
      Z<-Z[,colnames(Z)!="(Intercept)",drop=FALSE]
      colnames(Z)<-savecolnamesZ
    }
  } else {
    termsSD<-NULL
  }
  X<-X[Keep, ,drop=FALSE]
  Y<-Y[Keep]
  if (!is.null(weights)){
    weights<-weights[Keep]
    robust<-TRUE
    X<-X[!is.na(weights), ,drop=FALSE]
    Y<-Y[!is.na(weights)]
    if (!is.null(formulaSD) & !SameModelMEANSD){
      Z<-Z[!is.na(weights), ,drop=FALSE]
    }
    weights<-weights[!is.na(weights)]
    weighted<-TRUE
    weights<-weights*length(weights)/sum(weights)
  } else {
    weights<-rep_len(1,length(Y))
    weighted<-FALSE
  }
  NoVarModData<-data.frame(Y,weights)

  No.Obs<-length(Y)
  termsMODEL<-list(termsMEAN,termsSD)
  formulaMODEL<-list('meaneq' = formulaMEAN, 'sdeq' = formulaSD)
  # If specified that there is no constant in the model remove it from the data frame.
  if (!constantMEAN){
    savecolnames<-colnames(X)[colnames(X)!="(Intercept)"]
    X<-X[,colnames(X)!="(Intercept)",drop=FALSE]
    colnames(X)<-savecolnames
  }

  # collect variable means and check which variables are binary.
  XVarMeans<-apply(X,2,mean)
  XVarBinary<-apply(X,2,.checkbinary)

  Heteroskedastic<-TRUE # set model to heteroskedastic, check call and then switch if not.
  if (!is.null(formulaSD)){
    # will check if the formula for the mean is identical to the
    # the right hand side variables of a formula are accessible via terms in the attribute term.labels
    meaneqnames<-attr(terms(as.formula(formulaMEAN)),"term.labels")
    sdeqnames<-attr(terms(as.formula(formulaSD)),"term.labels")
    if (sum(is.na(match(meaneqnames,sdeqnames)))==sum(is.na(match(sdeqnames,meaneqnames))) & sum(is.na(match(sdeqnames,meaneqnames)))==0){
      if (constantSD==constantMEAN){
        SameModelMEANSD<-TRUE
        # collect the names and column numbers of variables that are in both the mean and variance equation
        meanandvarNAME<-colnames(X)
        meanandvarLOC<-c(1:ncol(X))
        meanandvarLOCZ<-meanandvarLOC<-meanandvarLOC[meanandvarNAME!="(Intercept)"]
        meanandvarNAME<-meanandvarNAME[meanandvarNAME!="(Intercept)"]
        BothMeanVar<-data.frame(meanandvarNAME,meanandvarLOC,meanandvarLOCZ,stringsAsFactors=FALSE)
        ZVarMeans<-XVarMeans
        ZVarBinary<-XVarBinary
      }
    }
  }

  if (!is.null(formulaSD) & !SameModelMEANSD){
    if (!constantSD){Z<-Z[,colnames(Z)!="(Intercept)",drop=FALSE]}
    ZVarMeans<-apply(Z,2,mean)
    ZVarBinary<-apply(Z,2,.checkbinary)
    # collect the names and column numbers of variables that are in both the mean and variance equation
    # find the colnames of Z that are the same as the colnames of X
    meanandvarLOC<-c(1:ncol(X))[!is.na(match(colnames(X),colnames(Z)))]
    # find the columns of Z that are in Z and X
    meanandvarLOCZ<-match(colnames(X),colnames(Z))[!is.na(match(colnames(X),colnames(Z)))]
    meanandvarNAME<-colnames(X)[meanandvarLOC]
    meanandvarLOC<-meanandvarLOC[meanandvarNAME!="(Intercept)"]
    meanandvarLOCZ<-meanandvarLOCZ[meanandvarNAME!="(Intercept)"]
    meanandvarNAME<-meanandvarNAME[meanandvarNAME!="(Intercept)"]
    BothMeanVar<-data.frame(meanandvarNAME,meanandvarLOC,meanandvarLOCZ,stringsAsFactors=FALSE)
  } else if (is.null(formulaSD) & !SameModelMEANSD){
    Z<-as.matrix(rep(1,nrow(X)),ncol=1)
    ZVarMeans<-NULL
    ZVarBinary<-NULL
    Heteroskedastic<-FALSE
    if (is.null(delta) & is.null(threshparam)){
      # if no formula for the standard deviation is given, the threshold parameters are not specified or the standard deviation is not specified then use the unit variance assumption.
      calcdelta<-function(x){eval({z<-x;sdmodel})-1}
      delta<-uniroot(calcdelta,c(-10,10),extendInt="yes",tol=.Machine$double.eps)$root # solve for delta to get the unit variance assumption
    }
    BothMeanVar=NULL
  }


  checkoutcomes<-.checkoutcomes(Y,Force=Force)
  listoutcomes<-checkoutcomes[[1]]
  no.outcomes<-checkoutcomes[[2]]
  #return(list(Y,X,X,weights,link,sdmodel,beta,delta,threshparam,analhessian,robust,start))
  output<-oglmx.fit(Y,X,Z,w=weights,
                    link = link,sdmodel = sdmodel,
                    beta=beta,delta=delta,
                    threshparam=threshparam,
                    analhessian=analhessian,
                    robustmatrix=robust,
                    start=start,
                    savemodelframe=savemodelframe,
                    optmeth = optmeth)
  output<-append(output,list(call=call,terms=termsMODEL,formula=formulaMODEL,NoVarModData=NoVarModData,Hetero=Heteroskedastic,BothEq=BothMeanVar,varMeans=list(XVarMeans,ZVarMeans),varBinary=list(XVarBinary,ZVarBinary)))
  class(output)<-"oglmx"
  output
}
