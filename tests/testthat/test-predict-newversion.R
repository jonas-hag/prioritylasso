if (!requireNamespace("rlang", quietly = TRUE)) {
  stop("The 'rlang' package is required for running this test. Please install it.")
}

# compare the current predict.prioritylasso version with the results from the
# CRAN version
# pl3 and pl3a are excluded because in the original prioritylasso function,
# an error was fixed, so the result is now different (use correct family for
# unpenalized first block and cvoffset)
set.seed(1234)
x_data <- matrix(rnorm(50*500),50,500)
# with single missing values
x_data_single_miss <- x_data
x_data_single_miss[1, 30] <- NA
# blockwise missing
x_data_miss_b1 <- x_data
x_data_miss_b1[1:20, 1:75] <- NA
x_data_miss_b2 <- x_data
x_data_miss_b2[21:30, 76:200] <- NA
x_data_miss_b12 <- x_data
x_data_miss_b12[30:40, 1:200] <- NA
x_data_miss_b123 <- x_data
x_data_miss_b123[30:40, 1:200] <- NA
x_data_miss_b123[41:45, 201:500] <- NA
x_data_comp_mis <- x_data
x_data_comp_mis[1, 1:500] <- NA
pl1a <- prioritylasso(X = x_data, Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:75, block2=76:200, block3=201:500),
                      block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5)
pred_pl1a <- predict(pl1a, newdata = x_data, type = "response")

set.seed(1234)
pl1b <- prioritylasso(X = x_data, Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:15, block2=16:200, block3=201:500),
                      block1.penalization = FALSE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5)
pred_pl1b <- predict(pl1b, newdata = x_data, type = "response")

set.seed(1234)
pl1_miss_comp <- suppressWarnings(prioritylasso(X = x_data_miss_b2, Y = rnorm(50), family = "gaussian", type.measure = "mse",
                                                blocks = list(block1=1:75, block2=76:200, block3=201:500),
                                                block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5,
                                                mcontrol = missing.control(handle.missingdata = "impute.offset")))

set.seed(1234)
pl1_miss_ign <- suppressWarnings(prioritylasso(X = x_data_miss_b12, Y = rnorm(50), family = "gaussian", type.measure = "mse",
                                               blocks = list(block1=1:75, block2=76:200, block3=201:500),
                                               block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5,
                                               mcontrol = missing.control(handle.missingdata = "ignore")))

set.seed(1234)
pl1_miss_avlb <- suppressWarnings(prioritylasso(X = x_data_miss_b123, Y = rnorm(50), family = "gaussian", type.measure = "mse",
                                                blocks = list(block1=1:75, block2=76:200, block3=201:500),
                                                block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5,
                                                mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                                           impute.offset.cases = "available.cases")))

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
pl3b <- suppressWarnings(prioritylasso(X = x_data, Y = rbinom(n=50, size=1, prob=0.5), family = "binomial",
                                       type.measure = "auc", blocks = list(block1=1:45, block2=46:200, block3=201:500),
                                       block1.penalization = TRUE, standardize = TRUE, nfolds = 4,
                                       cvoffset = TRUE))
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

################################################################################
# generate models that have information in the outcome
# gaussian outcome
set.seed(4973)
x_data_info <- cbind(matrix(rnorm(200 * 50), nrow = 200, ncol = 50),
                     matrix(rnorm(200 * 50), nrow = 200, ncol = 50),
                     matrix(rnorm(200 * 50), nrow = 200, ncol = 50))

x_data_info_miss_b1 <- x_data_info
x_data_info_miss_b1[1:20, 1:50] <- NA
# take the first variables (of each sd version) to have an influence on the
# outcome
y_data <- x_data_info[, c(1:2)] %*%
  matrix(c(10, -9), ncol = 1) +
  x_data_info[, 51:60] %*%
  matrix(c(1, -1, 2, -2, 2.5, 1.5, -1.5, 6, -2, 3), ncol = 1) +
  x_data_info[, c(101:105)] %*%
  matrix(c(0.5, 2, -4.8, 2, 7), ncol = 1)

y_data_gauss <- y_data +
  matrix(rnorm(200), ncol = 1)

# binomial outcome
prob_data <- 1 / (1 + exp(-y_data))
y_data_bin <- rbinom(200, 1, prob_data)

pl_gauss_info <- prioritylasso(X = x_data_info,
                               Y = y_data_gauss,
                               family = "gaussian",
                               type.measure = "mse",
                               blocks = list(block1 = 1:50, block2 = 51:100, block3 = 101:150),
                               block1.penalization = TRUE,
                               lambda.type = "lambda.1se",
                               standardize = TRUE,
                               nfolds = 5)

pl_gauss_info_miss_comp <- suppressWarnings(prioritylasso(X = x_data_info_miss_b1,
                                                          Y = y_data_gauss,
                                                          family = "gaussian",
                                                          type.measure = "mse",
                                                          blocks = list(block1 = 1:50, block2 = 51:100, block3 = 101:150),
                                                          block1.penalization = TRUE,
                                                          lambda.type = "lambda.1se",
                                                          standardize = TRUE,
                                                          nfolds = 5,
                                                          mcontrol = missing.control(handle.missingdata = "impute.offset")))

pl_binom_info <- prioritylasso(X = x_data_info,
                               Y = y_data_bin,
                               family = "binomial",
                               type.measure = "class",
                               blocks = list(block1 = 1:50, block2 = 51:100, block3 = 101:150),
                               block1.penalization = TRUE,
                               lambda.type = "lambda.1se",
                               standardize = TRUE,
                               nfolds = 5)


# load the comparison data from the CRAN version
path <- test_path("data_predict_test", "pred_pl1_cran.rds")
pred_pl1_cran <- readRDS(path)
path <- test_path("data_predict_test", "pred_pl1a_cran.rds")
pred_pl1a_cran <- readRDS(path)
path <- test_path("data_predict_test", "pred_pl1b_cran.rds")
pred_pl1b_cran <- readRDS(path)
path <- test_path("data_predict_test", "pred_pl2_cran.rds")
pred_pl2_cran <- readRDS(path)
path <- test_path("data_predict_test", "pred_pl2a_cran.rds")
pred_pl2a_cran <- readRDS(path)
path <- test_path("data_predict_test", "pred_pl2b_cran.rds")
pred_pl2b_cran <- readRDS(path)
path <- test_path("data_predict_test", "pred_pl3b_cran.rds")
pred_pl3b_cran <- readRDS(path)
path <- test_path("data_predict_test", "pred_pl4_cran.rds")
pred_pl4_cran <- readRDS(path)
path <- test_path("data_predict_test", "pred_pl5a_cran.rds")
pred_pl5a_cran <- readRDS(path)
path <- test_path("data_predict_test", "pred_pl5b_cran.rds")
pred_pl5b_cran <- readRDS(path)


context("tests for the further developed predict.prioritylasso")

test_that("the further developed predict.prioritylasso leads to the same results as the CRAN version", {
  
  expect_equal(pred_pl1, pred_pl1_cran)
  expect_equal(pred_pl1a, pred_pl1a_cran)
  expect_equal(pred_pl1b, pred_pl1b_cran)
  expect_equal(pred_pl2, pred_pl2_cran)
  expect_equal(pred_pl2a, pred_pl2a_cran)
  expect_equal(pred_pl2b, pred_pl2b_cran)
  expect_equal(pred_pl3b, pred_pl3b_cran)
  expect_equal(pred_pl4, pred_pl4_cran)
  expect_equal(pred_pl5a, pred_pl5a_cran)
  expect_equal(pred_pl5b, pred_pl5b_cran)
  
})

# write own expect function to check for differences
expect_difference <- function(object, comparison) {
  # 1. Capture object and label
  act <- quasi_label(rlang::enquo(object), arg = "object")
  act_2 <- quasi_label(rlang::enquo(comparison), arg = "object")
  
  # 2. Call expect()
  expect(
    !isTRUE(all.equal(act$val, act_2$val)),
    sprintf("%s is equal to %s.", act$lab, act_2$lab)
  )
  
  # 3. Invisibly return the value
  invisible(act$val)
}

# for NA
expect_na <- function(object) {
  # 1. Capture object and label
  act <- quasi_label(rlang::enquo(object), arg = "object")
  
  # 2. Call expect()
  expect(
    is.na(act$val) && length(act$val) == 1,
    sprintf("%s is not NA but of class %s or is not of length 1.",
            act$lab, class(act$val))
  )
  
  # 3. Invisibly return the value
  invisible(act$val)
  
}

# no NA
expect_no_na <- function(object) {
  # 1. Capture object and label
  act <- quasi_label(rlang::enquo(object), arg = "object")
  
  # 2. Call expect()
  expect(
    sum(is.na(as.vector(act$val))) == 0,
    sprintf("%s contains at least one NA.", act$lab)
  )
  
  # 3. Invisibly return the value
  invisible(act$val)
}

test_that("errors are thrown due to wrong input for object", {
  expect_error(predict(x_data, newdata = x_data, type = "response"))
})

pl_no_data <- prioritylasso(X = x_data, Y = rnorm(50), family = "gaussian",
                            type.measure = "mse",
                            blocks = list(block1=1:75, block2=76:200, block3=201:500),
                            block1.penalization = TRUE,
                            lambda.type = "lambda.1se", standardize = TRUE,
                            nfolds = 5, return.x = FALSE)

test_that("an error is thrown when no x value is available", {
  expect_error(predict(pl_no_data))
})

test_that("newdata is handled correctly", {
  # newdata can be a matrix or a data.frame
  expect_is(predict(pl1a, newdata = x_data, type = "response"), "matrix")
  expect_is(predict(pl1a, newdata = as.data.frame(x_data), type = "response"),
            "matrix")
  # newdata needn't to be supplied
  expect_is(predict(pl1a, type = "response"), "matrix")
  # can't have the wrong dimension
  expect_error(predict(pl1a, newdata = x_data[, -4], type = "response"))
  # can't be a vector
  expect_error(predict(pl1a, newdata = 1:500, type = "response"))
  # can't have a fully missing observation
  expect_error(predict(pl1a, newdata = x_data_comp_mis, type = "response"))
  expect_error(predict(pl1a, newdata = x_data_comp_mis, type = "response",
                       handle.missingtestdata = "omit.prediction"))
})

test_that("type is handled correctly", {
  # can only be link or response
  expect_is(predict(pl1a, newdata = x_data, type = "response"), "matrix")
  expect_is(predict(pl1a, newdata = x_data, type = "link"), "matrix")
  expect_error(predict(pl1a, newdata = x_data, type = "dumm"))
  # for gaussian response it should be the same
  expect_equal(predict(pl1a, newdata = x_data, type = "response"),
               predict(pl1a, newdata = x_data, type = "link"))
  # for binomial outcome it should be different
  expect_difference(predict(pl3b, newdata = x_data, type = "response"),
                    predict(pl3b, newdata = x_data, type = "link"))
})

test_that("handle.missingtestdata is handled correctly", {
  expect_error(predict(pl1a, newdata = x_data, type = "response",
                       handle.missingtestdata = "qual"))
  
  # single missing values are not accepted for any parameter value
  expect_error(predict(pl1a, newdata = x_data_single_miss, type = "response",
                       handle.missingtestdata = "none"))
  expect_error(predict(pl1a, newdata = x_data_single_miss, type = "response",
                       handle.missingtestdata = "omit.prediction"))
  expect_error(predict(pl1a, newdata = x_data_single_miss, type = "response",
                       handle.missingtestdata = "set.zero"))
  expect_error(predict(pl1a, newdata = x_data_single_miss, type = "response",
                       handle.missingtestdata = "impute.block"))
  
  # for complete data, you have to use "none"
  expect_error(predict(pl1a, newdata = x_data, type = "response",
                       handle.missingtestdata = "omit.prediction"))
  expect_error(predict(pl1a, newdata = x_data, type = "response",
                       handle.missingtestdata = "set.zero"))
  expect_error(predict(pl1a, newdata = x_data, type = "response",
                       handle.missingtestdata = "impute.block"))
  
  # for omit.prediction, you can use blockwise missing data
  expect_is(predict(pl1a, newdata = x_data_miss_b1, type = "response",
                    handle.missingtestdata = "omit.prediction"), "matrix")
  expect_is(predict(pl1a, newdata = x_data_miss_b12, type = "response",
                    handle.missingtestdata = "omit.prediction"), "matrix")
  expect_is(predict(pl1a, newdata = x_data_miss_b123, type = "response",
                    handle.missingtestdata = "omit.prediction"), "matrix")
  
  # if used, then this should lead to NA in the response
  expect_na(predict(pl1a, newdata = x_data_miss_b1, type = "response",
                    handle.missingtestdata = "omit.prediction")[1])
  
  # set.zero can be used also when trained on complete data
  expect_is(predict(pl1a, newdata = x_data_miss_b1, type = "response",
                    handle.missingtestdata = "set.zero"), "matrix")
  expect_is(predict(pl1a, newdata = x_data_miss_b12, type = "response",
                    handle.missingtestdata = "set.zero"), "matrix")
  expect_is(predict(pl1a, newdata = x_data_miss_b123, type = "response",
                    handle.missingtestdata = "set.zero"), "matrix")
  
  # set.zero can be used when trained on missing data, but only with missing
  # test data
  expect_is(predict(pl1_miss_ign, newdata = x_data_miss_b2, type = "response",
                    handle.missingtestdata = "set.zero"),
            "matrix")
  expect_no_na(predict(pl1_miss_ign, newdata = x_data_miss_b2, type = "response",
                       handle.missingtestdata = "set.zero"))
  
  # impute.block complete
  expect_is(predict(pl1_miss_comp, newdata = x_data_miss_b1, type = "response",
                    handle.missingtestdata = "impute.block"), "matrix")
  expect_no_na(predict(pl1_miss_comp, newdata = x_data_miss_b1, type = "response",
                       handle.missingtestdata = "impute.block"))
  
  # impute.block available
  expect_is(suppressWarnings(predict(pl1_miss_avlb, newdata = x_data_miss_b1, type = "response",
                                     handle.missingtestdata = "impute.block")), "matrix")
  expect_is(suppressWarnings(predict(pl1_miss_avlb, newdata = x_data_miss_b2, type = "response",
                                     handle.missingtestdata = "impute.block")), "matrix")
  expect_is(predict(pl1_miss_avlb, newdata = x_data_miss_b12, type = "response",
                    handle.missingtestdata = "impute.block"), "matrix")
  expect_is(predict(pl1_miss_avlb, newdata = x_data_miss_b123, type = "response",
                    handle.missingtestdata = "impute.block"), "matrix")
  expect_no_na(predict(pl1_miss_avlb, newdata = x_data_miss_b123, type = "response",
                       handle.missingtestdata = "impute.block"))
  
})

test_that("include.allintercepts is handled correctly", {
  # input check
  expect_error(predict(pl_gauss_info, newdata = x_data_info, type = "response",
                       include.allintercepts = "dumm"))
  
  # try out with gaussian and binomial outcome
  # use a model which acutally contains information for the outcome, not only
  # random outcome (otherwise the intercepts will be 0)
  expect_difference(predict(pl_gauss_info, newdata = x_data_info, type = "response"),
                    predict(pl_gauss_info, newdata = x_data_info, type = "response",
                            include.allintercepts = TRUE))
  expect_difference(predict(pl_binom_info, newdata = x_data_info, type = "response"),
                    predict(pl_binom_info, newdata = x_data_info, type = "response",
                            include.allintercepts = TRUE))
  
  # for survival outcome, it should be the same (because the cox model doesn't
  # have intercepts)
  expect_equal(predict(pl4, newdata = x, type = "response"),
               predict(pl4, newdata = x, type = "response",
                       include.allintercepts = TRUE))
  
  # try out with missing data
  expect_difference(predict(pl_gauss_info, newdata = x_data_info_miss_b1, type = "response",
                            handle.missingtestdata = "set.zero"),
                    predict(pl_gauss_info, newdata = x_data_info_miss_b1, type = "response",
                            handle.missingtestdata = "set.zero",
                            include.allintercepts = TRUE))
  expect_difference(predict(pl_gauss_info_miss_comp, newdata = x_data_info_miss_b1, type = "response",
                            handle.missingtestdata = "impute.block"),
                    predict(pl_gauss_info_miss_comp, newdata = x_data_info_miss_b1, type = "response",
                            handle.missingtestdata = "impute.block",
                            include.allintercepts = TRUE))
  expect_difference(predict(pl1_miss_avlb, newdata = x_data_miss_b123, type = "response",
                            handle.missingtestdata = "impute.block"),
                    predict(pl1_miss_avlb, newdata = x_data_miss_b123, type = "response",
                            handle.missingtestdata = "impute.block",
                            include.allintercepts = TRUE))
})

test_that("use.blocks is handled correctly", {
  # input checks
  expect_is(predict(pl1a, newdata = x_data, type = "response",
                    use.blocks = "all"), "matrix")
  expect_error(predict(pl1a, newdata = x_data, type = "response",
                       use.blocks = "dumm"))
  # no blocks other than are present in the data
  expect_error(predict(pl1a, newdata = x_data, type = "response",
                       use.blocks = 4))
  expect_error(predict(pl1a, newdata = x_data, type = "response",
                       use.blocks = -1))
  expect_error(predict(pl1a, newdata = x_data, type = "response",
                       use.blocks = 1:4))
  # there should be differences if only a subset of blocks are used
  expect_equal(predict(pl_gauss_info, newdata = x_data_info, type = "response",
                       use.blocks = "all"),
               predict(pl_gauss_info, newdata = x_data_info, type = "response",
                       use.blocks = 1:3))
  expect_difference(predict(pl_gauss_info, newdata = x_data_info, type = "response",
                            use.blocks = "all"),
                    predict(pl_gauss_info, newdata = x_data_info, type = "response",
                            use.blocks = 1))
  expect_difference(predict(pl_gauss_info, newdata = x_data_info, type = "response",
                            use.blocks = "all"),
                    predict(pl_gauss_info, newdata = x_data_info, type = "response",
                            use.blocks = 1:2))
  expect_difference(predict(pl_gauss_info, newdata = x_data_info, type = "response",
                            use.blocks = "all"),
                    predict(pl_gauss_info, newdata = x_data_info, type = "response",
                            use.blocks = c(1, 3)))
  
  # also works with missing data
  expect_difference(predict(pl_gauss_info_miss_comp, newdata = x_data_info_miss_b1, type = "response",
                            handle.missingtestdata = "impute.block"),
                    predict(pl_gauss_info_miss_comp, newdata = x_data_info_miss_b1, type = "response",
                            handle.missingtestdata = "impute.block",
                            use.blocks = c(1, 3)))
  
  # also works with use.allintercepts
  expect_difference(predict(pl_gauss_info, newdata = x_data_info, type = "response",
                            use.blocks = "all",
                            use.allintercepts = TRUE),
                    predict(pl_gauss_info, newdata = x_data_info, type = "response",
                            use.blocks = c(1, 3),
                            use.allintercepts = TRUE))
})