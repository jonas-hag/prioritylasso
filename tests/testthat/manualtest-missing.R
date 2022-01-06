set.seed(1234)
x_origin <- matrix(rnorm(50*500),50,500)
x_m1 <- x_origin
x_m1[1:10, 76:200] <- NA
pl_m1 <- prioritylasso(X = x_m1,
                      Y = rnorm(50),
                      family = "gaussian",
                      type.measure = "mse",
                      blocks = list(block1=1:75, block2=76:200, block3=201:500),
                      block1.penalization = TRUE,
                      lambda.type = "lambda.1se",
                      standardize = TRUE,
                      nfolds = 5,
                      # handle.missingdata = "ignore")
                      mcontrol = missing.control(handle.missingdata = "impute.offset"))

set.seed(1234)
pl_m1b <- prioritylasso(X = x_m1,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:75, block2=76:200, block3=201:500),
                       block1.penalization = TRUE,
                       lambda.type = "lambda.min",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

x_m2 <- x_origin
x_m2[1:10, 201:500] <- NA
set.seed(1234)
pl_m2 <- prioritylasso(X = x_m2,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:75, block2=76:200, block3=201:500),
                       block1.penalization = TRUE,
                       lambda.type = "lambda.1se",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

set.seed(1234)
pl_m2b<- prioritylasso(X = x_m2,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:75, block2=76:200, block3=201:500),
                       block1.penalization = TRUE,
                       lambda.type = "lambda.min",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

x_m3 <- x_origin
x_m3[1:10, 76:200] <- NA
x_m3[11:20, 201:500] <- NA
set.seed(1234)
pl_m3 <- prioritylasso(X = x_m3,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:75, block2=76:200, block3=201:500),
                       block1.penalization = TRUE,
                       lambda.type = "lambda.1se",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

set.seed(1234)
pl_m3b <- prioritylasso(X = x_m3,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:75, block2=76:200, block3=201:500),
                       block1.penalization = TRUE,
                       lambda.type = "lambda.min",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

x_m7 <- x_origin
x_m7[1:10, 1:75] <- NA
set.seed(1234)
pl_m7 <- prioritylasso(X = x_m7,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:75, block2=76:200, block3=201:500),
                       block1.penalization = TRUE,
                       lambda.type = "lambda.1se",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

set.seed(1234)
pl_m7b <- prioritylasso(X = x_m7,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:75, block2=76:200, block3=201:500),
                       block1.penalization = TRUE,
                       lambda.type = "lambda.min",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

x_m8 <- x_origin
x_m8[1:10, 1:75] <- NA
x_m8[11:20, 76:200] <- NA
x_m8[21:30, 201:500] <- NA
set.seed(1234)
pl_m8 <- prioritylasso(X = x_m8,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:75, block2=76:200, block3=201:500),
                       block1.penalization = TRUE,
                       lambda.type = "lambda.1se",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

set.seed(1234)
pl_m8b <- prioritylasso(X = x_m8,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:75, block2=76:200, block3=201:500),
                       block1.penalization = TRUE,
                       lambda.type = "lambda.min",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

x_m4 <- x_origin
x_m4[1:10, 16:200] <- NA
set.seed(1234)
pl_m4 <- prioritylasso(X = x_m4,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:15, block2=16:200, block3=201:500),
                       block1.penalization = FALSE,
                       lambda.type = "lambda.1se",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

set.seed(1234)
pl_m4b <- prioritylasso(X = x_m4,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:15, block2=16:200, block3=201:500),
                       block1.penalization = FALSE,
                       lambda.type = "lambda.min",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

x_m5 <- x_origin
x_m5[1:10, 201:500] <- NA
set.seed(1234)
pl_m5 <- prioritylasso(X = x_m5,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:15, block2=16:200, block3=201:500),
                       block1.penalization = FALSE,
                       lambda.type = "lambda.1se",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

set.seed(1234)
pl_m5b <- prioritylasso(X = x_m5,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:15, block2=16:200, block3=201:500),
                       block1.penalization = FALSE,
                       lambda.type = "lambda.min",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

x_m6 <- x_origin
x_m6[1:10, 16:200] <- NA
x_m6[11:20, 201:500] <- NA
set.seed(1234)
pl_m6 <- prioritylasso(X = x_m6,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:15, block2=16:200, block3=201:500),
                       block1.penalization = FALSE,
                       lambda.type = "lambda.1se",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

set.seed(1234)
pl_m6b <- prioritylasso(X = x_m6,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:15, block2=16:200, block3=201:500),
                       block1.penalization = FALSE,
                       lambda.type = "lambda.min",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")


x_m9 <- x_origin
x_m9[1:10, 1:15] <- NA
set.seed(1234)
pl_m9 <- prioritylasso(X = x_m9,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:15, block2=16:200, block3=201:500),
                       block1.penalization = FALSE,
                       lambda.type = "lambda.1se",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

set.seed(1234)
pl_m9b <- prioritylasso(X = x_m9,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:15, block2=16:200, block3=201:500),
                       block1.penalization = FALSE,
                       lambda.type = "lambda.min",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

x_m10 <- x_origin
x_m10[1:10, 1:15] <- NA
x_m10[11:20, 16:200] <- NA
x_m10[21:30, 201:500] <- NA
set.seed(1234)
pl_m10 <- prioritylasso(X = x_m10,
                       Y = rnorm(50),
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1=1:15, block2=16:200, block3=201:500),
                       block1.penalization = FALSE,
                       lambda.type = "lambda.1se",
                       standardize = TRUE,
                       nfolds = 5,
                       handle.missingdata = "ignore")

set.seed(1234)
pl_m10b <- prioritylasso(X = x_m10,
                        Y = rnorm(50),
                        family = "gaussian",
                        type.measure = "mse",
                        blocks = list(block1=1:15, block2=16:200, block3=201:500),
                        block1.penalization = FALSE,
                        lambda.type = "lambda.min",
                        standardize = TRUE,
                        nfolds = 5,
                        handle.missingdata = "ignore")

set.seed(4973)
# generate test data with different variances (and mean)
train_data <- cbind(matrix(rnorm(200 * 50), nrow = 200, ncol = 50),
                    matrix(rnorm(200 * 50), nrow = 200, ncol = 50),
                    matrix(rnorm(200 * 50), nrow = 200, ncol = 50))
# take the first variables (of each sd version) to have an influence on the
# outcome
train_y <- train_data[, c(1:2)] %*%
  matrix(c(10, -9), ncol = 1) +
  train_data[, 51:60] %*%
  matrix(c(1, -1, 2, -2, 2.5, 1.5, -1.5, 6, -2, 3), ncol = 1) +
  train_data[, c(101:105)] %*%
  matrix(c(0.5, 2, -4.8, 2, 7), ncol = 1) +
  matrix(rnorm(200), ncol = 1)

# delete some covariables
train_data[1:10, 1:50] <- NA
train_data[11:20, 51:100] <- NA
train_data[21:30, 101:150] <- NA

set.seed(1234)
pl_m11 <- prioritylasso(X = train_data,
                        Y = train_y,
                        family = "gaussian",
                        type.measure = "mse",
                        blocks = list(block1 = 1:50, block2 = 51:100, block3 = 101:150),
                        block1.penalization = TRUE,
                        lambda.type = "lambda.1se",
                        standardize = TRUE,
                        nfolds = 5,
                        mcontrol = missing.control(handle.missingdata = "impute.offset"))

# this should produce an error
train_data_2 <- train_data
train_data_2[11:20, 101:150] <- NA
pl_m11b <- prioritylasso(X = train_data_2,
                        Y = train_y,
                        family = "gaussian",
                        type.measure = "mse",
                        blocks = list(block1 = 1:50, block2 = 51:100, block3 = 101:150),
                        block1.penalization = TRUE,
                        lambda.type = "lambda.1se",
                        standardize = TRUE,
                        nfolds = 5,
                        mcontrol = missing.control(handle.missingdata = "impute.offset"))

################################################################################
# test use of 0 or intercept for handle.missingdata = ignore
pl_m12 <- prioritylasso(X = train_data,
                        Y = train_y,
                        family = "gaussian",
                        type.measure = "mse",
                        blocks = list(block1 = 1:50, block2 = 51:100, block3 = 101:150),
                        block1.penalization = TRUE,
                        lambda.type = "lambda.1se",
                        standardize = TRUE,
                        nfolds = 5,
                        mcontrol = missing.control(handle.missingdata = "ignore"))
pl_m12b <- prioritylasso(X = train_data,
                        Y = train_y,
                        family = "gaussian",
                        type.measure = "mse",
                        blocks = list(block1 = 1:50, block2 = 51:100, block3 = 101:150),
                        block1.penalization = FALSE,
                        lambda.type = "lambda.1se",
                        standardize = TRUE,
                        nfolds = 5,
                        mcontrol = missing.control(handle.missingdata = "ignore",
                                                   offset.firstblock = "intercept"))
pl_m12c <- prioritylasso(X = train_data,
                         Y = train_y,
                         family = "gaussian",
                         type.measure = "mse",
                         blocks = list(block1 = 1:50, block2 = 51:100, block3 = 101:150),
                         block1.penalization = TRUE,
                         lambda.type = "lambda.1se",
                         standardize = TRUE,
                         nfolds = 5,
                         mcontrol = missing.control(handle.missingdata = "ignore",
                                                    offset.firstblock = "intercept"))
