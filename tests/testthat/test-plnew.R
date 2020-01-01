# because for the unpenalised first block, when gaussian or binomial is used,
# the glm saves the data in an environment (it is not explicitly passed as a
# data frame), the environments differ between different objects, even when
# the same seed is used. This leads to fails in the tests, therefore the part
# block1unpen$data is set to NULL in these cases
# this information is not used by any function, so it can be omitted
set.seed(1234)
pl1a <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:75, block2=76:200, block3=201:500),
                      block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5)

set.seed(1234)
pl1b <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:15, block2=16:200, block3=201:500),
                      block1.penalization = FALSE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5)
pl1b$block1unpen$data <- NULL

###

set.seed(1234)
pl1 <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                     blocks = list(block1=1:75, block2=76:200, block3=201:500), max.coef = c(5,5,5),
                     block1.penalization = TRUE, lambda.type = "lambda.min", standardize = TRUE, nfolds = 5,
                     cvoffset = TRUE)

set.seed(1234)
pl2 <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                     blocks = list(block1=1:15, block2=16:200, block3=201:500),
                     block1.penalization = FALSE, lambda.type = "lambda.min", standardize = TRUE, nfolds = 5,
                     cvoffset = TRUE)
pl2$block1unpen$data <- NULL

set.seed(1234)
pl2a <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:20, block2=21:200, block3=201:500),
                      max.coef = c(Inf,Inf,Inf), block1.penalization = FALSE, lambda.type = "lambda.1se", nfolds = 5,
                      cvoffset = TRUE)
pl2a$block1unpen$data <- NULL

set.seed(1234)
pl2b <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:20, block2=21:200, block3=201:500),
                      max.coef = c(Inf,Inf,Inf), block1.penalization = TRUE, lambda.type = "lambda.1se", nfolds = 5,
                      cvoffset = TRUE)

set.seed(1234)
pl3 <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rbinom(n=50, size=1, prob=0.5), family = "binomial",
                     type.measure = "class", blocks = list(block1=1:15, block2=16:200, block3=201:500),
                     block1.penalization = FALSE, standardize = TRUE, nfolds = 5,
                     cvoffset = TRUE)
pl3$block1unpen$data <- NULL

set.seed(1234)
pl3a <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rbinom(n=50, size=1, prob=0.5), family = "binomial",
                     type.measure = "auc", blocks = list(block1=1:15, block2=16:200, block3=201:500),
                     block1.penalization = FALSE, standardize = TRUE, nfolds = 4,
                     cvoffset = TRUE)
pl3a$block1unpen$data <- NULL

set.seed(1234)
pl3b <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rbinom(n=50, size=1, prob=0.5), family = "binomial",
                      type.measure = "auc", blocks = list(block1=1:45, block2=46:200, block3=201:500),
                      block1.penalization = TRUE, standardize = TRUE, nfolds = 4,
                      cvoffset = TRUE)


### testing weights and foldid

set.seed(1234)
pl5a <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:75, block2=76:200, block3=201:500),
                      block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5,
                      weights = sample(rep(seq(1:2), length = 50)))

set.seed(1234)
pl5b <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:75, block2=76:200, block3=201:500),
                      block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 10,
                      foldid=sample(rep(seq(10),length=50)))

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

# load the comparison data from the CRAN version
if (grepl("prioritylasso$", getwd())) {
  path <- "data/cranversion_results/"
} else {
  path <- "../../data/cranversion_results/"
}

pl1_cran <- readRDS(paste0(path, "pl1_cran.rds"))
pl1a_cran <- readRDS(paste0(path, "pl1a_cran.rds"))
pl1b_cran <- readRDS(paste0(path, "pl1b_cran.rds"))
pl2_cran <- readRDS(paste0(path, "pl2_cran.rds"))
pl2a_cran <- readRDS(paste0(path, "pl2a_cran.rds"))
pl2b_cran <- readRDS(paste0(path, "pl2b_cran.rds"))
pl3_cran <- readRDS(paste0(path, "pl3_cran.rds"))
pl3a_cran <- readRDS(paste0(path, "pl3a_cran.rds"))
pl3b_cran <- readRDS(paste0(path, "pl3b_cran.rds"))
pl4_cran <- readRDS(paste0(path, "pl4_cran.rds"))
pl5a_cran <- readRDS(paste0(path, "pl5a_cran.rds"))
pl5b_cran <- readRDS(paste0(path, "pl5b_cran.rds"))

pl1b_cran$block1unpen$data <- NULL
pl2_cran$block1unpen$data <- NULL
pl2a_cran$block1unpen$data <- NULL
pl3_cran$block1unpen$data <- NULL
pl3a_cran$block1unpen$data <- NULL


library(testthat)

context("tests for the further developed prioritylasso")

test_that("the further developed prioritylasso leads to the same results as the CRAN version", {

  expect_equal(pl1, pl1_cran)
  expect_equal(pl1a, pl1a_cran)
  expect_equal(pl1b, pl1b_cran)
  expect_equal(pl2, pl2_cran)
  expect_equal(pl2a, pl2a_cran)
  expect_equal(pl2b, pl2b_cran)
  expect_equal(pl3, pl3_cran)
  expect_equal(pl3a, pl3a_cran)
  expect_equal(pl3b, pl3b_cran)
  expect_equal(pl4, pl4_cran)
  expect_equal(pl5a, pl5a_cran)
  expect_equal(pl5b, pl5b_cran)

})
