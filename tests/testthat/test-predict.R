set.seed(1234)
pl_bin1 <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rbinom(50,1,0.5), family = "binomial",
                       type.measure = "auc", blocks = list(block1=1:13,block2=14:200, block3=201:500),
                       block1.penalization = TRUE, lambda.type = "lambda.min", standardize = FALSE,
                       nfolds = 5)

newdata_bin1 <- matrix(rnorm(30*500),30,500)

bin1_predict <- predict(object = pl_bin1, newdata = newdata_bin1, type = "response")

###

set.seed(1234)
pl_bin2 <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rbinom(50,1,0.5), family = "binomial",
                         type.measure = "auc", blocks = list(block1=1:13,block2=14:200, block3=201:500),
                         block1.penalization = FALSE, lambda.type = "lambda.min", standardize = FALSE,
                         nfolds = 5)

newdata_bin2 <- matrix(rnorm(30*500),30,500)

bin2_predict <- predict(object = pl_bin2, newdata = newdata_bin2, type = "response")

###

n <- 50;p <- 300
nzc <- trunc(p/10)
set.seed(1234)
x <- matrix(rnorm(n*p), n, p)
beta <- rnorm(nzc)
fx <- x[, seq(nzc)]%*%beta/3
hx <- exp(fx)
# survival times:
ty <- rexp(n,hx)
# censoring indicator:
tcens <- rbinom(n = n,prob = .3,size = 1)
library(survival)
y <- Surv(ty, 1-tcens)
blocks <- list(bp1=1:20, bp2=21:200, bp3=201:300)
# run prioritylasso:
set.seed(1234)
pl_surv1 <- prioritylasso(x, y, family = "cox", type.measure = "deviance", blocks = blocks,
              block1.penalization = TRUE, lambda.type = "lambda.min", standardize = TRUE,
              nfolds = 5)

newdata_surv1 <- matrix(rnorm(30*300), 30, 300)

cox1_predict <- predict(object = pl_surv1, newdata = newdata_surv1, type = "response")

###

set.seed(1234)
pl_surv2 <- prioritylasso(x, y, family = "cox", type.measure = "deviance", blocks = blocks,
                          block1.penalization = FALSE, lambda.type = "lambda.min", standardize = TRUE,
                          nfolds = 5)

newdata_surv2 <- matrix(rnorm(30*300), 30, 300)

cox2_predict <- predict(object = pl_surv2, newdata = newdata_surv2, type = "response")

###

set.seed(1234)
pl_gaussian1 <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian",
                             type.measure = "mse", blocks = list(bp1=1:20, bp2=21:200, bp3=201:500),
                             max.coef = c(Inf,8,5), block1.penalization = TRUE,
                             lambda.type = "lambda.min", standardize = TRUE, nfolds = 5, cvoffset = TRUE)

newdata_gaussian1 <- matrix(rnorm(30*500),30,500)
newdata_gaussian2 <- matrix(rnorm(30*500),30,500)
newdata_gaussian2[1:10, 21:200] <- NA

set.seed(1234)
data_missing <- matrix(rnorm(50*500),50,500)
data_missing[1:10, 21:200] <- NA
pl_gaussian2 <- suppressWarnings(prioritylasso(X = data_missing, Y = rnorm(50), family = "gaussian",
                              type.measure = "mse", blocks = list(bp1=1:20, bp2=21:200, bp3=201:500),
                              max.coef = c(Inf,8,5), block1.penalization = TRUE,
                              lambda.type = "lambda.min", standardize = TRUE, nfolds = 5, cvoffset = FALSE,
                              mcontrol = missing.control(handle.missingdata = "ignore")))



gaussian1_predict <- predict(object = pl_gaussian1, newdata = newdata_gaussian1, type = "link")


###


library(testthat)

context("tests for predict.prioritylasso")

test_that("testing predictions for binomial family", {

  expect_true(dim(bin1_predict)[1] == 30)
  expect_true(class(bin2_predict)[1] == "matrix")

})

test_that("testing predictions for gaussian family", {

  expect_true(dim(gaussian1_predict)[2] == 1)

})

test_that("also data with missing entries is used for prediction", {
  expect_true(dim(predict(pl_gaussian1, newdata = newdata_gaussian2,
                          handle.missingtestdata = "set.zero"))[1] == 30)
})

test_that("prediction works with object trained with missing data", {
  expect_true(dim(predict(pl_gaussian2, newdata = newdata_gaussian1))[1] == 30)
})



test_that("testing predictions for cox family", {

  expect_true(dim(cox1_predict)[1] == 30)
  expect_true(dim(cox2_predict)[2] == 1)

})

test_that("testing other stuff", {
  expect_that(predict(object = pl_gaussian1, newdata = newdata_gaussian1, type = ""),
              throws_error("'arg' should be one of \"link\", \"response\""))
})
