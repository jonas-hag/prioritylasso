# compare the current predict.prioritylasso version with the results from the
# CRAN version
set.seed(1234)
x_data <- matrix(rnorm(50*500),50,500)
pl1a <- prioritylasso(X = x_data, Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:75, block2=76:200, block3=201:500),
                      block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5)
pred_pl1a <- predict(pl1a, newdata = x_data, type = "response")

set.seed(1234)
pl1b <- prioritylasso(X = x_data, Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:15, block2=16:200, block3=201:500),
                      block1.penalization = FALSE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5)
pred_pl1b <- predict(pl1b, newdata = x_data, type = "response")

###

set.seed(1234)
pl1 <- prioritylasso(X = x_data, Y = rnorm(50), family = "gaussian", type.measure = "mse",
                     blocks = list(block1=1:75, block2=76:200, block3=201:500), max.coef = c(5,5,5),
                     block1.penalization = TRUE, lambda.type = "lambda.min", standardize = TRUE, nfolds = 5,
                     cvoffset = TRUE)
pred_pl1 <- predict(pl1, newdata = x_data, type = "response")

set.seed(1234)
pl2 <- prioritylasso(X = x_data, Y = rnorm(50), family = "gaussian", type.measure = "mse",
                     blocks = list(block1=1:15, block2=16:200, block3=201:500),
                     block1.penalization = FALSE, lambda.type = "lambda.min", standardize = TRUE, nfolds = 5,
                     cvoffset = TRUE)
pred_pl2 <- predict(pl2, newdata = x_data, type = "response")

set.seed(1234)
pl2a <- prioritylasso(X = x_data, Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:20, block2=21:200, block3=201:500),
                      max.coef = c(Inf,Inf,Inf), block1.penalization = FALSE, lambda.type = "lambda.1se", nfolds = 5,
                      cvoffset = TRUE)
pred_pl2a <- predict(pl2a, newdata = x_data, type = "response")

set.seed(1234)
pl2b <- prioritylasso(X = x_data, Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:20, block2=21:200, block3=201:500),
                      max.coef = c(Inf,Inf,Inf), block1.penalization = TRUE, lambda.type = "lambda.1se", nfolds = 5,
                      cvoffset = TRUE)
pred_pl2b <- predict(pl2b, newdata = x_data, type = "response")

set.seed(1234)
pl3 <- prioritylasso(X = x_data, Y = rbinom(n=50, size=1, prob=0.5), family = "binomial",
                     type.measure = "class", blocks = list(block1=1:15, block2=16:200, block3=201:500),
                     block1.penalization = FALSE, standardize = TRUE, nfolds = 5,
                     cvoffset = TRUE)
pred_pl3 <- predict(pl3, newdata = x_data, type = "response")

set.seed(1234)
pl3a <- prioritylasso(X = x_data, Y = rbinom(n=50, size=1, prob=0.5), family = "binomial",
                      type.measure = "auc", blocks = list(block1=1:15, block2=16:200, block3=201:500),
                      block1.penalization = FALSE, standardize = TRUE, nfolds = 4,
                      cvoffset = TRUE)
pred_pl3a <- predict(pl3a, newdata = x_data, type = "response")

set.seed(1234)
pl3b <- prioritylasso(X = x_data, Y = rbinom(n=50, size=1, prob=0.5), family = "binomial",
                      type.measure = "auc", blocks = list(block1=1:45, block2=46:200, block3=201:500),
                      block1.penalization = TRUE, standardize = TRUE, nfolds = 4,
                      cvoffset = TRUE)
pred_pl3b <- predict(pl3b, newdata = x_data, type = "response")


### testing weights and foldid

set.seed(1234)
pl5a <- prioritylasso(X = x_data, Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:75, block2=76:200, block3=201:500),
                      block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5,
                      weights = sample(rep(seq(1:2), length = 50)))
pred_pl5a <- predict(pl5a, newdata = x_data, type = "response")

set.seed(1234)
pl5b <- prioritylasso(X = x_data, Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:75, block2=76:200, block3=201:500),
                      block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 10,
                      foldid=sample(rep(seq(10),length=50)))
pred_pl5b <- predict(pl5b, newdata = x_data, type = "response")

###

# cox
n <- 50;p <- 300
nzc <- trunc(p/10)
set.seed(1234)
x <- matrix(rnorm(n*p), n, p)
beta <- rnorm(nzc)
fx <- x[, seq(nzc)]%*%beta/3
hx <- exp(fx)
ty <- rexp(n,hx)
tcens <- rbinom(n = n,prob = .3,size = 1)
library(survival)
y <- Surv(ty, 1-tcens)
blocks <- list(block1=1:20, block2=21:200, block3=201:300)

set.seed(1234)
pl4 <- prioritylasso(x, y, family = "cox", type.measure = "deviance", blocks = blocks,
                     block1.penalization = FALSE, lambda.type = "lambda.min", standardize = TRUE, nfolds = 5,
                     cvoffset = TRUE)
pred_pl4 <- predict(pl4, newdata = x, type = "response")

# load the comparison data from the CRAN version
  path <- "../../data/data_predict_test/"


pred_pl1_cran <- readRDS(paste0(path, "pred_pl1_cran.rds"))
pred_pl1a_cran <- readRDS(paste0(path, "pred_pl1a_cran.rds"))
pred_pl1b_cran <- readRDS(paste0(path, "pred_pl1b_cran.rds"))
pred_pl2_cran <- readRDS(paste0(path, "pred_pl2_cran.rds"))
pred_pl2a_cran <- readRDS(paste0(path, "pred_pl2a_cran.rds"))
pred_pl2b_cran <- readRDS(paste0(path, "pred_pl2b_cran.rds"))
pred_pl3_cran <- readRDS(paste0(path, "pred_pl3_cran.rds"))
pred_pl3a_cran <- readRDS(paste0(path, "pred_pl3a_cran.rds"))
pred_pl3b_cran <- readRDS(paste0(path, "pred_pl3b_cran.rds"))
pred_pl4_cran <- readRDS(paste0(path, "pred_pl4_cran.rds"))
pred_pl5a_cran <- readRDS(paste0(path, "pred_pl5a_cran.rds"))
pred_pl5b_cran <- readRDS(paste0(path, "pred_pl5b_cran.rds"))

library(testthat)

context("tests for the further developed predict.prioritylasso")

test_that("the further developed predict.prioritylasso leads to the same results as the CRAN version", {
  
  expect_equal(pred_pl1, pred_pl1_cran)
  expect_equal(pred_pl1a, pred_pl1a_cran)
  expect_equal(pred_pl1b, pred_pl1b_cran)
  expect_equal(pred_pl2, pred_pl2_cran)
  expect_equal(pred_pl2a, pred_pl2a_cran)
  expect_equal(pred_pl2b, pred_pl2b_cran)
  expect_equal(pred_pl3, pred_pl3_cran)
  expect_equal(pred_pl3a, pred_pl3a_cran)
  expect_equal(pred_pl3b, pred_pl3b_cran)
  expect_equal(pred_pl4, pred_pl4_cran)
  expect_equal(pred_pl5a, pred_pl5a_cran)
  expect_equal(pred_pl5b, pred_pl5b_cran)
  
})