testthat::context("oglm implements ordered and interval regression models")


# I - CREATE DATA -------------------------

set.seed(242)
n<-250
x1<-sample(c(0,1),n,replace=TRUE,prob=c(0.75,0.25))
x2<-vector("numeric",n)
x2[x1==0]<-sample(c(0,1),n-sum(x1==1),replace=TRUE,prob=c(2/3,1/3))
z<-rnorm(n,0.5)
# create latent outcome variable
latenty<-0.5+1.5*x1-0.5*x2+0.5*z+rnorm(n,sd=exp(0.5*x1-0.5*x2))
# observed y has four possible values: -1,0,1,2
# threshold values are: -0.5, 0.5, 1.5.
y<-vector("numeric",n)
y[latenty< -0.5]<--1
y[latenty>= -0.5 & latenty<0.5]<- 0
y[latenty>= 0.5 & latenty<1.5]<- 1
y[latenty>= 1.5]<- 2
dataset<-data.frame(y,x1,x2)


requireNamespace("oglmx", quietly = TRUE)




# POSSIBILITY TO CHANGE OPTIMIZATION METHOD -------------


testthat::test_that("Default method to Newton-Raphson", {
  mod1 <- oglm::oglmx(y ~ x1 + x2 + z,
                      data=dataset,link="probit",constantMEAN=FALSE,
                      constantSD=FALSE,delta=0,threshparam=NULL)
  mod2 <- oglm::oglmx(y ~ x1 + x2 + z,
                      data=dataset,link="probit",constantMEAN=FALSE,
                      constantSD=FALSE,delta=0,threshparam=NULL,
                      optmeth = "NR")
  mod1$call <- NULL
  mod2$call <- NULL

  testthat::expect_identical(
    mod1, mod2
  )

})


testthat::test_that("Raise error when argument not matched", {

  testthat::expect_error(
    oglm::oglmx(y ~ x1 + x2 + z,
                data=dataset,link="probit",constantMEAN=FALSE,
                constantSD=FALSE,delta=0,threshparam=NULL,
                optmeth = "nr")
  )
  testthat::expect_error(
    oglm::oglmx(y ~ x1 + x2 + z,
                data=dataset,link="probit",constantMEAN=FALSE,
                constantSD=FALSE,delta=0,threshparam=NULL,
                optmeth = "anything")
  )


})



# II - TEST OPTIONS -------------------

# II-1- FORMULA =====

testthat::test_that(
  "Excepting call and formula, elements are the same", {

    formulastr <- oglm::oglmx("y ~ x1 + x2",
                              data=dataset,link="probit",constantMEAN=FALSE,
                              constantSD=FALSE,delta=0,threshparam=NULL)
    formulaformal <- oglm::oglmx(y ~ x1 + x2,
                                 data=dataset,link="probit",constantMEAN=FALSE,
                                 constantSD=FALSE,delta=0,threshparam=NULL)

    formulastr$formula <- NULL
    formulastr$call <- NULL
    formulaformal$formula <- NULL
    formulaformal$call <- NULL

    testthat::expect_equal(
      formulastr,
      formulaformal
    )

  })


testthat::test_that(
  "Same call produces an error with oglmx package", {
    oglm::oglmx("y ~ x1 + x2",
                data=dataset,link="probit",constantMEAN=FALSE,
                constantSD=FALSE,delta=0,threshparam=NULL)
    testthat::expect_error(
      oglmx::oglmx("y ~ x1 + x2",
                   data=dataset,link="probit",constantMEAN=FALSE,
                   constantSD=FALSE,delta=0,threshparam=NULL)
    )
})



# COMPARE TO STATA ---------------------


requireNamespace("haven", quietly = TRUE)

stata <- haven::read_dta("http://www.stata-press.com/data/r13/womenwage.dta")

wage_lbounds <- sort(as.numeric(
  unique(as.character(stata$wagecat)))
)
wage_lbounds <- wage_lbounds[-length(wage_lbounds)]


model <- oglm::oglmx(
  data = stata,
  formula = "wagecat ~ age + I(age^2) + nev_mar + rural + school + tenure",
  formulaSD = NULL,
  link = "probit",
  constantMEAN = TRUE,
  constantSD = TRUE,
  threshparam =  log(wage_lbounds)
)
# see Stata rintreg manual p.7


coefs_stata <- c(.7084023, .0645589, -.0010812, -.0058151, -.2098361,
                 .0804832, .0397144, -.906989)
# log(sigma) also estimated
sd_stata <- c(.3593193, .0249954, .0004115, .0454867,
              .0439454, .0076783, .0058001,
              .0356265)

testthat::expect_equal(
  coefs_stata,
  as.numeric(model$coefficients),
  tolerance = 1e-5
)


