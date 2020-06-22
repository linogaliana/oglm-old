testthat::context("gumbel distribution well implemented")


testthat::test_that("rgumbel and Rgumbel generate random values", {

  testthat::expect_equal(
    length(rgumbel(100L)),
    100L
  )

  testthat::expect_equal(
    length(rGumbel(100L)),
    100L
  )

})

x <- rGumbel(100L)
mu <- 1L
beta <- 2L

testthat::test_that("pgumbel and pGumbel", {

  testthat::expect_true(
    min(pgumbel(x, loc = mu, scale = beta)) >= 0L
  )

  testthat::expect_true(
    max(pgumbel(x, loc = mu, scale = beta)) <= 1L
  )

  testthat::expect_true(
    min(pGumbel(x, loc = mu, scale = beta)) >= 0L
  )

  testthat::expect_true(
    max(pGumbel(x, loc = mu, scale = beta)) <= 1L
  )



  testthat::expect_true(
    min(pgumbel(x, loc = mu, scale = beta, lower.tail = FALSE)) >= 0L
  )

  testthat::expect_true(
    max(pgumbel(x, loc = mu, scale = beta, lower.tail = FALSE)) <= 1L
  )

  testthat::expect_true(
    min(pGumbel(x, loc = mu, scale = beta, lower.tail = FALSE)) >= 0L
  )

  testthat::expect_true(
    max(pGumbel(x, loc = mu, scale = beta, lower.tail = FALSE)) <= 1L
  )


})




testthat::test_that("dgumbel and dGumbel", {

  testthat::expect_equal(
    sum(dgumbel(x, loc = mu, scale = beta) >= rep(0L, length(x))),
    length(x)
  )

  testthat::expect_equal(
    sum(dGumbel(x, loc = mu, scale = beta) >= rep(0L, length(x))),
    length(x)
  )

  testthat::expect_true(
    max(dgumbel(x, loc = mu, scale = beta, log = TRUE))<=0L
  )

  testthat::expect_true(
    max(dGumbel(x, loc = mu, scale = beta, log = TRUE))<=0L
  )

})
