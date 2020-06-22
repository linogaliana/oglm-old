#' Predictions from ordered discrete model
#'
#' @inheritParams stats::predict.glm
#' @param object An ordered discrete model computed using `oglmx` package
#' @param type Prediction considered. Can be \emph{'class'}
#'  (predicted label as maximum probability label, default);
#'  \emph{'probs'} (probabilities of each observation to belong to a particular class)
#'  or \emph{'latent'} (\eqn{Y^*} latent variable)
#' @return Depends of the value of prediction \code{type} \describe{
#'   \item{class}{Most likely label, i.e. \eqn{l = \arg \max_l p_j}}
#'   \item{probs}{Prediced probabilities for each label, i.e. \eqn{p_ij} matrix with
#'    $i$ observation index and $j$ class label}
#'   \item{latent}{Prediced value in latent space \eqn{y^*}}
#' }
#'
#'
#' @importFrom stats .checkMFClasses delete.response
#' @importFrom stats predict
#' @importFrom stats model.frame model.matrix napredict pcauchy
#' @importFrom stats plogis pnorm rcauchy rlogis rnorm terms
#' @export


predict.oglmx <- function (object, newdata = NULL, type = c("class", "probs","latent","xb"), ...){

  # CHECK IF oglmx OBJECT
  # --------------------------
  if (!inherits(object, "oglmx")) stop("not a \"oglmx\" object")

  # CHECK PREDICTION TYPE
  # --------------------------
  type <- match.arg(type)

  # FIT IF missing(newdata)
  # ---------------------------

  if (missing(newdata)){
    #Y <- object$fitted
    stop('fit method not yet implemented: use predict')
  }

  newdata <- as.data.frame(newdata)
  formula <- object$formula$`meaneq`


  # CREATE TERMS TO REPLICATE MASS::predict.polr BEHAVIOR
  # --------------------------------------------------

  # Transform formula in terms
  object$terms <- terms(formula)

  # Keep only covariates
  Terms <- delete.response(object$terms)

  # Covariates matrix
  m <- model.frame(Terms, newdata, na.action = function(x) x,
                   xlev = object$factorvars)

  # Check factors are ok
  if (!is.null(cl <- attr(Terms, "dataClasses")))
    .checkMFClasses(cl, m)

  X <- model.matrix(Terms, m, contrasts = object$contrasts)

  xnocol <- match(names(object$coefficients),
                  colnames(X),
                  nomatch = 0L)
  if (length(xnocol)>0L){
    X2 <- X[, xnocol]
  }


  # Finalize covariates matrix by removing intercept
  xint <- match("(Intercept)", colnames(X), nomatch = 0L)
  if (xint > 0L){
    X <- X[, -xint, drop = FALSE]
  }

  # Parameters for logistic/normal/... distribution
  n <- nrow(X)
  q <- length(object$allparams$threshparam)


  # MAKE PREDICTION
  # -------------------------------

  coeff_list <- object$coefficients

  # REMOVE THRESHOLD COEFFICIENTS
  if (sum(grepl("Threshold", names(coeff_list)))>0){
    coeff <- coeff_list[-grep("Threshold", names(coeff_list))]
  }else{
    coeff <- coeff_list
  }

  # REMOVE STANDARD ERROR RELATED PARAMETER
  coeff <- coeff[!is.na(names(coeff))]

  # POSSIBILITY THAT SOME FIXED EFFECTS ARE NOT PRESENT IN X2
  # REMOVE THEM FROM coeff
  coeff2 <- coeff[names(coeff) %in% colnames(X2)]

  # x*\beta vector (nb: intercept column should be added back)
  eta <- drop(X2 %*% coeff2)

  if (type == "xb") return(eta)

  # if type == "latent", we simulate epsilon
  if (type == "latent"){

    epsilon_distribution <- switch(object$link, logit = rlogis, probit = rnorm,
                                   loglog = rgumbel, cloglog = rGumbel, cauchit = rcauchy)


    # sigma identified to 1 in case of non-user defined interval regression
    if (sum(grepl("Threshold", names(coeff_list)))>0){

      simulated_epsilon <- epsilon_distribution(length(eta))
      sigma_vector <- rep(1L, times = length(eta)) # A vérifier

    } else{

      sigma_vector <- sigma(object, newdata = newdata)

      if (object$link == "logit"){
        m <- 0
        sigma_vector <- sqrt(3)*sigma_vector/pi #Variance for logistic is s^2 * pi^2/3
        simulated_epsilon <- sapply(sigma_vector, function(s) epsilon_distribution(1L, location = 0,
                                                                                   scale = s)
        )
      } else if (object$link == "probit"){
        simulated_epsilon <- sapply(sigma_vector, function(s) epsilon_distribution(1L, mean = 0,
                                                                                   sd = s)
        )
      } else{
        # Other distributions for epsilon
        simulated_epsilon <- epsilon_distribution(length(eta))
      }

    }

    eta <- eta + simulated_epsilon

  }

  # Which distribution should be applied ?
  pfun <- switch(object$link, logit = plogis, probit = pnorm,
                 loglog = pgumbel, cloglog = pGumbel, cauchit = pcauchy)

  # Transform from latent space y = x*beta to probabilities
  cumpr <- matrix(pfun(matrix(object$allparams$threshparam, n, q, byrow = TRUE) -
                         eta), , q)

  # CREATE PREDICTION
  # ---------------------------------

  Y <- t(apply(cumpr, 1L, function(x) diff(c(0, x, 1))))
  dimnames(Y) <- list(rownames(X), object$Outcomes)


  if (missing(newdata) && !is.null(object$na.action)){
    Y <- napredict(object$na.action, Y)
  }


  if (type == "class"){
    factor(max.col(Y), levels = seq_along(object$Outcomes),
           labels = object$Outcomes)
  }else if (type == "probs") {
    drop(Y)
  } else{
    # When type = "latent"
    return(list("xb" = eta - simulated_epsilon, "y_latent_pred" = eta,
                "sigma_vector" = sigma_vector))
  }
}
