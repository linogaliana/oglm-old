testthat::context("Formalize the construction of a variance model")


iris$y <- sample(1:5, size = nrow(iris),
                 replace = TRUE)

# I - MODEL WHERE THRESHOLDS ARE ESTIMATED ------

# A/ NO formulaSD ARGUMENT -----

ordered_logit <- REtage::ordered_model_threshold(
  iris,
  formula = "y ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width",
  link = "logit",
  thresholds = NULL, constantSD = FALSE,
  constantMEAN = FALSE, delta = 0)



testthat::test_that(
  "[no newdata] Model for variance is a column of 1s when formulaSD is NULL",{
    testthat::expect_equal(
      variance_model(ordered_logit), matrix(rep(1, nrow(iris)))
    )
  }
)

testthat::test_that(
  "[with newdata] Model for variance is a column of 1s when formulaSD is NULL",{
    testthat::expect_equal(
      variance_model(ordered_logit, newdata = data.frame(runif(100L))),
      matrix(rep(1, 100L))
    )
  }
)


# B/ MODEL HAS A formulaSD ARGUMENT

ordered_logit <- REtage::ordered_model_threshold(
  iris,
  formula = "y ~ Sepal.Length + Sepal.Width",
  formulaSD = as.formula("y ~ Sepal.Width"),
  link = "logit",
  thresholds = NULL,
  constantMEAN = FALSE, delta = 0)

testthat::expect_error(
  Z_no_newdata <- variance_model(ordered_logit)
)




# MODEL WHERE THRESHOLDS ARE KNOWN ------
