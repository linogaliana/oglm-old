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


# POSSIBILITY TO CHANGE OPTIMIZATION METHOD -------------


testthat::test_that("Default method to Newton-Raphson", {
  mod1 <- oglm::oglmx(y ~ x1 + x2 + z,
                      data=dataset,link="probit",constantMEAN=FALSE,
                      constantSD=FALSE,delta=0,threshparam=NULL)
  mod2 <- oglm::oglmx(y ~ x1 + x2 + z,
                      data=dataset,link="probit",constantMEAN=FALSE,
                      constantSD=FALSE,delta=0,threshparam=NULL,
                      optmeth = "NR")
  mod3 <- oglm::oglmx(y ~ x1 + x2 + z,
                      data=dataset,link="probit",constantMEAN=FALSE,
                      constantSD=FALSE,delta=0,threshparam=NULL,
                      optmeth = "nr")
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
