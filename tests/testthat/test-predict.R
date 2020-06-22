context("Predictions are equivalent to MASS::predict method for polr")

iris$y <- sample(1:5, size = nrow(iris),
                 replace = TRUE)



ordered_logit <- REtage::ordered_model_threshold(
  iris,
  formula = "y ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width",
  link = "logit",
  thresholds = NULL, constantSD = FALSE,
  constantMEAN = FALSE, delta = 0)

ordered_logit_MASS <- MASS::polr("I(factor(y)) ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width",
                                 data = iris,
                                 Hess = TRUE,
                                 method = "logistic")

ordered_probit <- REtage::ordered_model_threshold(
  iris,
  formula = "y ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width",
  link = "probit",
  thresholds = NULL, constantSD = FALSE,
  constantMEAN = FALSE, delta = 0)

ordered_probit_MASS <- MASS::polr("I(factor(y)) ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width",
                                  data = iris,
                                  Hess = TRUE,
                                  method = "probit")


# ERROR WHEN NEWDATA IS NULL -----------------------

testthat::expect_error(predict(ordered_probit), regexp = "fit method not yet implemented")
testthat::expect_error(predict(ordered_logit), regexp = "fit method not yet implemented")


# TESTS FOR PROBIT ------------------------

oglm_class <- predict(ordered_probit, iris)
polr_class <- predict(ordered_probit_MASS, iris)

testthat::test_that("Class predictions are the same with oglmx and polr",{
  testthat::expect_equal(
    as.numeric(as.character(polr_class)),
    as.numeric(as.character(oglm_class))
  )
})


oglm_probs <- predict(ordered_probit, iris, type = "probs")
polr_probs <- predict(ordered_probit_MASS, iris, type = "probs")

head(oglm_probs)
head(polr_probs)


testthat::test_that("Prediction matrix has expected dimensions",{
  testthat::expect_equal(
    dim(oglm_probs),
    dim(polr_probs)
  )
})


testthat::test_that("Prediction matrix has values close from MASS::polr [tolerance: at 1e-5 level]",{
  testthat::expect_equal(
    oglm_probs,
    polr_probs,
    tolerance = 2e-5
  )
})


# TESTS FOR LOGIT -------------------------

oglm_class <- predict(ordered_logit, iris)
polr_class <- predict(ordered_logit_MASS, iris)


testthat::test_that("Class predictions are the same with oglmx and polr",{
  testthat::expect_equal(as.numeric(as.character(polr_class)),
                         as.numeric(as.character(oglm_class))
  )
})


testthat::test_that("Prediction matrix has expected dimensions",{
  testthat::expect_equal(
    dim(oglm_probs),
    dim(polr_probs)
  )
})


testthat::test_that("Prediction matrix has values close from MASS::polr [tolerance: at 1e-5 level]",{
  testthat::expect_equal(
    oglm_probs,
    polr_probs,
    tolerance = 2e-5
  )
})



# TEST LATENT MODEL PREDICTION -----------------


## OUTPUT FORMAT AS EXPECTED

xb_logit <- predict(ordered_logit, iris, type = "latent")
xb_probit <- predict(ordered_probit, iris, type = "latent")



testthat::test_that("Return two elements in list: xb, y_latent_pred and sigma",{
  testthat::expect_equal(typeof(xb_logit), "list")
  testthat::expect_equal(typeof(xb_probit), "list")
  testthat::expect_equal(names(xb_logit), c("xb","y_latent_pred","sigma_vector"))
  testthat::expect_equal(names(xb_probit), c("xb","y_latent_pred","sigma_vector"))
}
)

testthat::test_that("Return three elements in list: xb, y_latent_pred and sigma",{
  testthat::expect_true(is.numeric(xb_logit$xb))
  testthat::expect_true(is.numeric(xb_probit$xb))
  testthat::expect_true(is.numeric(xb_probit$y_latent_pred))
  testthat::expect_true(is.numeric(xb_probit$y_latent_pred))
  #  expect_true(is.numeric(xb_logit$sigma))
  #  expect_true(is.numeric(xb_probit$sigma))
})


## CONSISTENCY BETWEEN LATENT AND CLASS MODEL PREDICTIONS --------------

set.seed(242)

n <- 250

x1 <- sample(c(0, 1), n,
             replace = TRUE,
             prob = c(0.75, 0.25))

x2 <- vector("numeric", n)
x2[x1 == 0] <- sample(c(0, 1),
                      n - sum(x1 == 1),
                      replace = TRUE,
                      prob = c(2/3, 1/3))

z <- rnorm(n, 0.5)
# create latent outcome variable
latenty <- -0.5 + 1.5 * x1 - 0.5 * x2 + 0.5 * z + rnorm(n, sd = exp(0.5 *
                                                                      x1 - 0.5 * x2))
# observed y has four possible values: -1,0,1,2
# threshold values are: -0.5, 0.5, 1.5.
y <- vector("numeric", n)
y[latenty < 0.5] <- -1
y[latenty >= 0.5 & latenty < 1.5] <- 0
y[latenty >= 1.5 & latenty < 2.5] <- 1
y[latenty >= 2.5] <- 2
dataset <- data.frame(y, x1, x2)

bounds <- c(0.5, 1.5, 2.5)
lbounds <- log(bounds)




ordered_probit <- REtage::ordered_model_threshold(
  dataset,
  formula = "y ~ x1+x2",
  link = "probit",
  thresholds = lbounds,
  constantSD = TRUE)


xb_probit <- predict(ordered_probit, dataset, type = "xb")
latent_probit <- predict(ordered_probit, newdata = dataset, type = "latent")




## MONTE-CARLO EXPERIMENT

# df <- data.frame(x = rnorm(10e5))
# df$y_latent   <- df$x + rnorm(nrow(df))
# df$y_observed <- Hmisc::cut2(df$y_latent,
#                                cuts = quantile(df$y_latent, probs = c(0.25,.5,.75)))
#
#
# ordered_probit <- ordered_model_threshold(
#   df,
#   formula = "y_observed ~ 0 + x",
#   link = "probit",
#   thresholds = quantile(df$y_latent, probs = c(0.25,.5,.75))
# )
#
# pred <- predict(ordered_probit, df, type = "latent")
# y_predict <- pred$y_latent_pred
# xb <- pred$xb
