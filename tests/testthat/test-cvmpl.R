# gaussian
cvm_pl1 <- cvm_prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                             blocks.list = list(list(block1=1:75,block2=76:200, block3=201:500), list(1:75, 201:500, 76:200)),
                             block1.penalization = TRUE, lambda.type = "lambda.min", standardize = TRUE, nfolds = 5)


cvm_pl1a <- cvm_prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                              blocks.list = list(list(block1=1:75,block2=76:200, block3=201:500), list(1:75, 201:500, 76:200)),
                              max.coef.list = list(c(10,5,7), c(10,7,3)),
                              block1.penalization = TRUE, lambda.type = "lambda.min", standardize = TRUE, nfolds = 5)


# binomial
cvm_pl2 <- cvm_prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rbinom(n=50, size=1, prob=0.5), family = "binomial",
                             type.measure = "auc", blocks.list = list(list(block1=1:75,block2=76:200, block3=201:500),
                                                                      list(1:75, 201:500, 76:200)),
                             block1.penalization = TRUE, lambda.type = "lambda.min",
                             standardize = TRUE, nfolds = 5)



# cox
n <- 50;p <- 300
nzc <- trunc(p/10)
x <- matrix(rnorm(n*p), n, p)
beta <- rnorm(nzc)
fx <- x[, seq(nzc)]%*%beta/3
hx <- exp(fx)
ty <- rexp(n,hx)
tcens <- rbinom(n = n,prob = .3,size = 1)
library(survival)
y <- Surv(ty, 1-tcens)
blocks <- list(list(block1=1:20, block2=21:200, block3=201:300), list(1:20, 201:300, 21:200))

cvm_pl3 <- cvm_prioritylasso(x, y, family = "cox", type.measure = "deviance", blocks.list = blocks,
                             block1.penalization = FALSE,
                             lambda.type = "lambda.min", standardize = TRUE, nfolds = 5)

# missingness
# -> checks that further arguments can be passed from cvm_prioritylasso to
# prioritylasso
set.seed(1234)
x_data <- matrix(rnorm(50*500),50,500)
x_data_miss_b123 <- x_data
x_data_miss_b123[30:40, 1:200] <- NA
x_data_miss_b123[41:45, 201:500] <- NA
set.seed(1234)
cvm_miss <- suppressWarnings(cvm_prioritylasso(
  X = x_data_miss_b123, Y = rnorm(50),
  family = "gaussian",
  type.measure = "mse",
  blocks.list = list(list(block1=1:75, block2=76:200, block3=201:500),
                     list(block1=76:200, block2=1:75, block3=201:500),
                     list(block1=201:500, block2=1:75, block3=76:200)),
  block1.penalization = TRUE,
  lambda.type = "lambda.1se",
  standardize = TRUE,
  nfolds = 5,
  mcontrol = missing.control(handle.missingdata = "impute.offset",
                             impute.offset.cases = "available.cases")))

test_that("testing error messages", {
  
  expect_that(cvm_prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                                blocks.list = list(list(1:75, 76:200, 201:500),list(1:75, 201:500)),
                                block1.penalization = TRUE, lambda.type = "lambda.min", standardize = TRUE, nfolds = 5),
              throws_error("Each predictor should be included in exactly one block."))
  
  expect_that(cvm_prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                                blocks.list = list(list(block1=1:75,block2=76:200, block3=201:500), list(1:75, 201:500, 76:200)),
                                max.coef.list = list(c(5,6,7)),
                                block1.penalization = TRUE, lambda.type = "lambda.min", standardize = TRUE, nfolds = 5),
              throws_error("blocks.list and max.coef.list must have the same length."))
  
  
  expect_that(cvm_prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian", type.measure = "mse",
                                blocks = list(list(block1=1:75,block2=76:200, block3=201:500), list(1:75, 76:500)),
                                max.coef.list = list(c(5,6,7), c(8,9,10)),
                                block1.penalization = TRUE, lambda.type = "lambda.min", standardize = TRUE, nfolds = 5),
              throws_error("blocks.list and the entries of max.coef.list must have the same length."))
  
  
})

test_that("testing cvm_prioritylasso", {
  
  expect_that(length(cvm_pl1$best.blocks), testthat::equals(3))
  expect_match(cvm_pl2$name, "AUC")
  expect_that(cvm_pl3, is_a("cvm_prioritylasso"))
  expect_true(cvm_pl1a$nzero[[1]] <= 10)
  expect_true(sum(unlist(cvm_pl1a$nzero)) <= 22)
  
})

test_that("cvm_prioritylasso contains the correct return elements", {
  expect_equal(names(cvm_miss), c("lambda.ind", "lambda.type", "lambda.min",
                                  "min.cvm", "nzero", "glmnet.fit",
                                  "name", "block1unpen", "best.blocks",
                                  "best.blocks.indices",
                                  "best.max.coef", "best.model","coefficients",
                                  "call"))
  expect_s3_class(cvm_miss[["best.model"]], "prioritylasso")
  expect_type(cvm_miss[["best.blocks.indices"]], "list")
})