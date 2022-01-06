# tests are done with browser/debugging in predict function

# check that all intercepts are correctly used
# unpenalised first block
set.seed(1234)
x_data <- matrix(rnorm(50*500),50,500)
pl1 <- prioritylasso(X = x_data, Y = rnorm(50), family = "gaussian", type.measure = "mse",
                      blocks = list(block1=1:15, block2=16:200, block3=201:500),
                      block1.penalization = FALSE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5)
pred_pl1 <- predict(pl1, newdata = x_data, type = "response",
                    include.allintercepts = TRUE)

# penalised first block
set.seed(1234)
x_data <- matrix(rnorm(50*500),50,500)
pl2 <- prioritylasso(X = x_data, Y = rnorm(50), family = "gaussian", type.measure = "mse",
                     blocks = list(block1=1:15, block2=16:200, block3=201:500),
                     block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5)
pred_pl2 <- predict(pl2, newdata = x_data, type = "response",
                    include.allintercepts = TRUE)

# check that y-backscaling is done
set.seed(1234)
x_data <- matrix(rnorm(50*500),50,500)
pl3 <- prioritylasso(X = x_data, Y = rnorm(50, mean = 3, sd = 2), family = "gaussian", type.measure = "mse",
                     blocks = list(block1=1:15, block2=16:200, block3=201:500),
                     block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5,
                     scale.y = TRUE)
pred_pl3 <- predict(pl3, newdata = x_data, type = "response",
                    include.allintercepts = TRUE)

################################################################################
# check that the imputation of missing test values works
set.seed(1234)
train_data <- cbind(matrix(rnorm(200 * 50, mean = 2, sd = 1), nrow = 200, ncol = 50),
                    matrix(rnorm(200 * 50, mean = 0, sd = 2), nrow = 200, ncol = 50),
                    matrix(rnorm(200 * 50, mean = -1, sd = 1), nrow = 200, ncol = 50))
# take the first variables (of each sd version) to have an influence on the
# outcome
train_y <- train_data[, c(1:2)] %*%
  matrix(c(10, -9), ncol = 1) +
  train_data[, 51:60] %*%
  matrix(c(1, -1, 2, -2, 2.5, 1.5, -1.5, 6, -2, 3), ncol = 1) +
  train_data[, 101:110] %*%
  matrix(c(0.5, 2, 4, -3, 1.5, 7, 0.7, 2, -1, 2.4), ncol = 1) +
  matrix(rnorm(200), ncol = 1)

set.seed(1234)
pl4 <- prioritylasso(X = train_data, Y = train_y, family = "gaussian", type.measure = "mse",
                     blocks = list(block1=1:50, block2=51:100, block3=101:150),
                     block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5,
                     mcontrol = missing.control(handle.missingdata = "impute.offset"))

test_data <- train_data
pred_pl4 <- predict(pl4, newdata = test_data, type = "response",
                    handle.missingtestdata = "none")
# should lead to an error
pred_pl4b <- predict(pl4, newdata = test_data, type = "response",
                    handle.missingtestdata = "impute.block")

test_data_miss <- test_data
test_data_miss[1:5, 1:50] <- NA
test_data_miss[6:10, 51:100] <- NA
test_data_miss[11:15, 101:150] <- NA

pred_pl4c <- predict(pl4, newdata = test_data_miss, type = "response",
                     handle.missingtestdata = "impute.block")

################################################################################
# use the data with missing entries in the train phase
set.seed(1234)
pl5 <- prioritylasso(X = test_data_miss, Y = train_y, family = "gaussian", type.measure = "mse",
                     blocks = list(block1=1:50, block2=51:100, block3=101:150),
                     block1.penalization = TRUE, lambda.type = "lambda.1se", standardize = TRUE, nfolds = 5,
                     mcontrol = missing.control(handle.missingdata = "impute.offset"))

pred_pl5 <- predict(pl5, newdata = test_data_miss, type = "response",
                    handle.missingtestdata = "impute.block")
