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
                      handle.missingdata = "ignore")

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