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

# delete some values
train_data_2 <- train_data
# in the first block
train_data_2[1:20, 1:50] <- NA
# in the second block
train_data_2[100:140, 51:100] <- NA
# in the third block
train_data_2[120:150, 101:150] <- NA


set.seed(1234)
pl_m1 <- prioritylasso(X = train_data_2,
                        Y = train_y,
                        family = "gaussian",
                        type.measure = "mse",
                        blocks = list(block1 = 1:50, block2 = 51:100, block3 = 101:150),
                        block1.penalization = TRUE,
                        lambda.type = "lambda.min",
                        standardize = TRUE,
                        nfolds = 5,
                        mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                   impute.offset.cases = "available.cases"))

pl_m2 <- prioritylasso(X = train_data,
                       Y = train_y,
                       family = "gaussian",
                       type.measure = "mse",
                       blocks = list(block1 = 1:50, block2 = 51:100, block3 = 101:150),
                       block1.penalization = TRUE,
                       lambda.type = "lambda.min",
                       standardize = TRUE,
                       nfolds = 5,
                       mcontrol = missing.control(handle.missingdata = "impute.offset",
                                                  impute.offset.cases = "available.cases"))

