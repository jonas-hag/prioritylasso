# because for the unpenalised first block, when gaussian or binomial is used,
# the glm saves the data in an environment (it is not explicitly passed as a
# data frame), the environments differ between different objects, even when
# the same seed is used. This leads to fails in the tests, therefore the part
# block1unpen$data is set to NULL in these cases
# this information is not used by any function, so it can be omitted
# also, all the names and dimnames of some components of the block1unpen are set
# to NULL because due to changed code, the names are now different (but contain
# the same values)
# also, the call/formula/terms are set to NULL in block1unpen because it changed
# due to different code (but still should lead to the same coefficients, which
# are tested)
# compared to the CRAN version, X and missing.data have been added in the value,
# delete it for comparison

set.seed(1234)
pl1a <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:75, block2=76:200, block3=201:500),
                      block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5)
pl1a$X <- NULL
pl1a$missing.data <- NULL

set.seed(1234)
pl1b <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:15, block2=16:200, block3=201:500),
                      block1.penalization = FALSE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5)

pl1b$X <- NULL
pl1b$missing.data <- NULL
pl1b$block1unpen$data <- NULL
names(pl1b$block1unpen$effects) <- NULL
dimnames(pl1b$block1unpen$R) <- NULL
dimnames(pl1b$block1unpen$qr$qr) <- NULL
names(pl1b$block1unpen$model) <- NULL
attr(pl1b$block1unpen$model, "terms") <- NULL
pl1b$block1unpen$call <- NULL
pl1b$block1unpen$formula <- NULL
pl1b$block1unpen$terms <- NULL

###

set.seed(1234)
pl1 <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                     blocks = list(block1=1:75, block2=76:200, block3=201:500), max.coef = c(5,5,5),
                     block1.penalization = TRUE, lambda.type = "lambda.min", standardize = TRUE, nfolds = 5,
                     cvoffset = TRUE)
pl1$X <- NULL
pl1$missing.data <- NULL

set.seed(1234)
pl2 <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                     blocks = list(block1=1:15, block2=16:200, block3=201:500),
                     block1.penalization = FALSE, lambda.type = "lambda.min", standardize = TRUE, nfolds = 5,
                     cvoffset = TRUE)
pl2$X <- NULL
pl2$missing.data <- NULL
pl2$block1unpen$data <- NULL
names(pl2$block1unpen$effects) <- NULL
dimnames(pl2$block1unpen$R) <- NULL
dimnames(pl2$block1unpen$qr$qr) <- NULL
names(pl2$block1unpen$model) <- NULL
attr(pl2$block1unpen$model, "terms") <- NULL
pl2$block1unpen$call <- NULL
pl2$block1unpen$formula <- NULL
pl2$block1unpen$terms <- NULL

set.seed(1234)
pl2a <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:20, block2=21:200, block3=201:500),
                      max.coef = c(Inf,Inf,Inf), block1.penalization = FALSE, lambda.type = "lambda.1se", nfolds = 5,
                      cvoffset = TRUE)
pl2a$X <- NULL
pl2a$missing.data <- NULL
pl2a$block1unpen$data <- NULL
names(pl2a$block1unpen$effects) <- NULL
dimnames(pl2a$block1unpen$R) <- NULL
dimnames(pl2a$block1unpen$qr$qr) <- NULL
names(pl2a$block1unpen$model) <- NULL
attr(pl2a$block1unpen$model, "terms") <- NULL
pl2a$block1unpen$call <- NULL
pl2a$block1unpen$formula <- NULL
pl2a$block1unpen$terms <- NULL

set.seed(1234)
pl2b <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:20, block2=21:200, block3=201:500),
                      max.coef = c(Inf,Inf,Inf), block1.penalization = TRUE, lambda.type = "lambda.1se", nfolds = 5,
                      cvoffset = TRUE)
pl2b$X <- NULL
pl2b$missing.data <- NULL

set.seed(1234)
pl3 <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rbinom(n=50, size=1, prob=0.5), family = "binomial",
                     type.measure = "class", blocks = list(block1=1:15, block2=16:200, block3=201:500),
                     block1.penalization = FALSE, standardize = TRUE, nfolds = 5,
                     cvoffset = TRUE)
pl3$X <- NULL
pl3$missing.data <- NULL
pl3$block1unpen$data <- NULL
names(pl3$block1unpen$effects) <- NULL
dimnames(pl3$block1unpen$R) <- NULL
dimnames(pl3$block1unpen$qr$qr) <- NULL
names(pl3$block1unpen$model) <- NULL
attr(pl3$block1unpen$model, "terms") <- NULL
pl3$block1unpen$call <- NULL
pl3$block1unpen$formula <- NULL
pl3$block1unpen$terms <- NULL

set.seed(1234)
pl3a <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rbinom(n=50, size=1, prob=0.5), family = "binomial",
                     type.measure = "auc", blocks = list(block1=1:15, block2=16:200, block3=201:500),
                     block1.penalization = FALSE, standardize = TRUE, nfolds = 4,
                     cvoffset = TRUE)
pl3a$X <- NULL
pl3a$missing.data <- NULL
pl3a$block1unpen$data <- NULL
names(pl3a$block1unpen$effects) <- NULL
dimnames(pl3a$block1unpen$R) <- NULL
dimnames(pl3a$block1unpen$qr$qr) <- NULL
names(pl3a$block1unpen$model) <- NULL
attr(pl3a$block1unpen$model, "terms") <- NULL
pl3a$block1unpen$call <- NULL
pl3a$block1unpen$formula <- NULL
pl3a$block1unpen$terms <- NULL

set.seed(1234)
pl3b <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rbinom(n=50, size=1, prob=0.5), family = "binomial",
                      type.measure = "auc", blocks = list(block1=1:45, block2=46:200, block3=201:500),
                      block1.penalization = TRUE, standardize = TRUE, nfolds = 4,
                      cvoffset = TRUE)
pl3b$X <- NULL
pl3b$missing.data <- NULL


### testing weights and foldid

set.seed(1234)
pl5a <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:75, block2=76:200, block3=201:500),
                      block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5,
                      weights = sample(rep(seq(1:2), length = 50)))
pl5a$X <- NULL
pl5a$missing.data <- NULL

set.seed(1234)
pl5b <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:75, block2=76:200, block3=201:500),
                      block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 10,
                      foldid=sample(rep(seq(10),length=50)))
pl5b$X <- NULL
pl5b$missing.data <- NULL

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
pl4$X <- NULL
pl4$missing.data <- NULL

names(pl4$block1unpen$means) <- NULL
pl4$block1unpen$terms <- NULL
names(pl4$block1unpen$assign) <- NULL
names(pl4$block1unpen$model) <- NULL
attr(pl4$block1unpen$model, "terms") <- NULL
pl4$block1unpen$formula <- NULL
pl4$block1unpen$call <- NULL



# load the comparison data from the CRAN version
path <- system.file("refactoring_results/pl1_cran.rds", package = "prioritylasso")
pl1_cran <- readRDS(path)
path <- system.file("refactoring_results/pl1a_cran.rds", package = "prioritylasso")
pl1a_cran <- readRDS(path)
path <- system.file("refactoring_results/pl1b_cran.rds", package = "prioritylasso")
pl1b_cran <- readRDS(path)
path <- system.file("refactoring_results/pl2_cran.rds", package = "prioritylasso")
pl2_cran <- readRDS(path)
path <- system.file("refactoring_results/pl2a_cran.rds", package = "prioritylasso")
pl2a_cran <- readRDS(path)
path <- system.file("refactoring_results/pl2b_cran.rds", package = "prioritylasso")
pl2b_cran <- readRDS(path)
path <- system.file("refactoring_results/pl3_cran.rds", package = "prioritylasso")
pl3_cran <- readRDS(path)
path <- system.file("refactoring_results/pl3a_cran.rds", package = "prioritylasso")
pl3a_cran <- readRDS(path)
path <- system.file("refactoring_results/pl3b_cran.rds", package = "prioritylasso")
pl3b_cran <- readRDS(path)
path <- system.file("refactoring_results/pl4_cran.rds", package = "prioritylasso")
pl4_cran <- readRDS(path)
path <- system.file("refactoring_results/pl5a_cran.rds", package = "prioritylasso")
pl5a_cran <- readRDS(path)
path <- system.file("refactoring_results/pl5b_cran.rds", package = "prioritylasso")
pl5b_cran <- readRDS(path)

pl1b_cran$block1unpen$data <- NULL
pl2_cran$block1unpen$data <- NULL
pl2a_cran$block1unpen$data <- NULL
pl3_cran$block1unpen$data <- NULL
pl3a_cran$block1unpen$data <- NULL

names(pl1b_cran$block1unpen$effects) <- NULL
dimnames(pl1b_cran$block1unpen$R) <- NULL
dimnames(pl1b_cran$block1unpen$qr$qr) <- NULL
names(pl1b_cran$block1unpen$model) <- NULL
attr(pl1b_cran$block1unpen$model, "terms") <- NULL
pl1b_cran$block1unpen$call <- NULL
pl1b_cran$block1unpen$formula <- NULL
pl1b_cran$block1unpen$terms <- NULL

names(pl2_cran$block1unpen$effects) <- NULL
dimnames(pl2_cran$block1unpen$R) <- NULL
dimnames(pl2_cran$block1unpen$qr$qr) <- NULL
names(pl2_cran$block1unpen$model) <- NULL
attr(pl2_cran$block1unpen$model, "terms") <- NULL
pl2_cran$block1unpen$call <- NULL
pl2_cran$block1unpen$formula <- NULL
pl2_cran$block1unpen$terms <- NULL

names(pl2a_cran$block1unpen$effects) <- NULL
dimnames(pl2a_cran$block1unpen$R) <- NULL
dimnames(pl2a_cran$block1unpen$qr$qr) <- NULL
names(pl2a_cran$block1unpen$model) <- NULL
attr(pl2a_cran$block1unpen$model, "terms") <- NULL
pl2a_cran$block1unpen$call <- NULL
pl2a_cran$block1unpen$formula <- NULL
pl2a_cran$block1unpen$terms <- NULL

names(pl3_cran$block1unpen$effects) <- NULL
dimnames(pl3_cran$block1unpen$R) <- NULL
dimnames(pl3_cran$block1unpen$qr$qr) <- NULL
names(pl3_cran$block1unpen$model) <- NULL
attr(pl3_cran$block1unpen$model, "terms") <- NULL
pl3_cran$block1unpen$call <- NULL
pl3_cran$block1unpen$formula <- NULL
pl3_cran$block1unpen$terms <- NULL

names(pl3a_cran$block1unpen$effects) <- NULL
dimnames(pl3a_cran$block1unpen$R) <- NULL
dimnames(pl3a_cran$block1unpen$qr$qr) <- NULL
names(pl3a_cran$block1unpen$model) <- NULL
attr(pl3a_cran$block1unpen$model, "terms") <- NULL
pl3a_cran$block1unpen$call <- NULL
pl3a_cran$block1unpen$formula <- NULL
pl3a_cran$block1unpen$terms <- NULL

names(pl4_cran$block1unpen$means) <- NULL
pl4_cran$block1unpen$terms <- NULL
names(pl4_cran$block1unpen$assign) <- NULL
names(pl4_cran$block1unpen$model) <- NULL
attr(pl4_cran$block1unpen$model, "terms") <- NULL
pl4_cran$block1unpen$formula <- NULL
pl4_cran$block1unpen$call <- NULL

# the calls to glmnet where changed (indexing of observations), but in these
# cases they should not change the calculations, therefore they are set to NULL

for (i in 1:3) {
  pl1_cran$glmnet.fit[[i]][["call"]] <- NULL
  pl1a_cran$glmnet.fit[[i]][["call"]] <- NULL
  pl2b_cran$glmnet.fit[[i]][["call"]] <- NULL
  pl3b_cran$glmnet.fit[[i]][["call"]] <- NULL
  
  pl5a_cran$glmnet.fit[[i]][["call"]] <- NULL
  pl5b_cran$glmnet.fit[[i]][["call"]] <- NULL
  
  pl1$glmnet.fit[[i]][["call"]] <- NULL
  pl1a$glmnet.fit[[i]][["call"]] <- NULL
  pl2b$glmnet.fit[[i]][["call"]] <- NULL
  pl3b$glmnet.fit[[i]][["call"]] <- NULL
  pl5a$glmnet.fit[[i]][["call"]] <- NULL
  pl5b$glmnet.fit[[i]][["call"]] <- NULL
}

# for the examples where block 1 is unpenalised, there are only 2 glmnet objects
for (i in 2:3) {
  pl1b_cran$glmnet.fit[[i]][["call"]] <- NULL
  pl2_cran$glmnet.fit[[i]][["call"]] <- NULL
  pl2a_cran$glmnet.fit[[i]][["call"]] <- NULL
  pl3_cran$glmnet.fit[[i]][["call"]] <- NULL
  pl3a_cran$glmnet.fit[[i]][["call"]] <- NULL
  pl4_cran$glmnet.fit[[i]][["call"]] <- NULL
  
  pl1b$glmnet.fit[[i]][["call"]] <- NULL
  pl2$glmnet.fit[[i]][["call"]] <- NULL
  pl2a$glmnet.fit[[i]][["call"]] <- NULL
  pl3$glmnet.fit[[i]][["call"]] <- NULL
  pl3a$glmnet.fit[[i]][["call"]] <- NULL
  pl4$glmnet.fit[[i]][["call"]] <- NULL
}

library(testthat)

context("tests for the further developed prioritylasso")

test_that("the further developed prioritylasso leads to the same results as the refactored version", {

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
