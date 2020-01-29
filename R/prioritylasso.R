#' Patient outcome prediction based on multi-omics data taking practitioners' preferences into account
#'
#' Fits successive Lasso models for several ordered blocks of (omics) data and takes the predicted values as an offset for the next block.
#'
#' For \code{block1.penalization = TRUE}, the function fits a Lasso model for each block. First, a standard Lasso for the first entry of \code{blocks} (block of priority 1) is fitted.
#' The predictions are then taken as an offset in the Lasso fit of the block of priority 2, etc.
#' For \code{block1.penalization = FALSE}, the function fits a model without penalty to the block of priority 1 (recommended as a block with clinical predictors where \code{p < n}).
#' This is either a generalized linear model for family "gaussian" or "binomial", or a Cox model. The predicted values are then taken as an offset in the following Lasso fit of the block with priority 2, etc. \cr
#'
#' The first entry of \code{blocks} contains the indices of variables of the block with priority 1 (first block included in the model).
#' Assume that \code{blocks = list(1:100, 101:200, 201:300)} then the block with priority 1 consists of the first 100 variables of the data matrix.
#' Analogously, the block with priority 2 consists of the variables 101 to 200 and the block with priority 3 of the variables 201 to 300.
#' 
#' \code{standardize = TRUE} leads to a standardisation of the covariables (\code{X}) in \code{glmnet} which is recommend by \code{glmnet}.
#' In case of an unpenalized first block, the covariables for the first block are not standardized.
#' Please note that the returned coefficients are rescaled to the original scale of the covariates as provided in \code{X}.
#' Therefore, new data in \code{predict.prioritylasso} should be on the same scale as \code{X}.
#' 
#' To use the method with blockwise missing data, one can set \code{handle.missingdata = ignore}.
#' Then, to calculate the coefficients for a given block only the observations with values for this blocks are used.
#' For the observations with missing values, the result from the previous block is used as the offset for the next block.
#' Crossvalidated offsets are not supported with \code{handle.missingdata = ignore}.
#' Please note that dealing with single missing values is not supported.
#' If only one value in a block for one observation is missing, still the complete observation for this block is treated as missing.
#' Normally, every observation gets a unique foldid which stays the same across all blocks for the call to \code{cv.glmnet}.
#' However when \code{handle.missingdata != none}, the foldid is set new for every block.
#'
#' @param X a (nxp) matrix of predictors with observations in rows and predictors in columns.
#' @param Y n-vector giving the value of the response (either continuous, numeric-binary 0/1, or \code{Surv} object).
#' @param weights observation weights. Default is 1 for each observation.
#' @param family should be "gaussian" for continuous \code{Y}, "binomial" for binary \code{Y}, "cox" for \code{Y} of type \code{Surv}.
#' @param type.measure accuracy/error measure computed in cross-validation. It should be "class" (classification error) or "auc" (area under the ROC curve) if \code{family="binomial"}, "mse" (mean squared error) if \code{family="gaussian"} and "deviance" if \code{family="cox"} which uses the partial-likelihood.
#' @param blocks list of the format \code{list(bp1=...,bp2=...,)}, where the dots should be replaced by the indices of the predictors included in this block. The blocks should form a partition of 1:p.
#' @param max.coef vector with integer values which specify the number of maximal coefficients for each block. The first entry is omitted if \code{block1.penalization = FALSE}. Default is \code{NULL}.
#' @param block1.penalization whether the first block should be penalized. Default is TRUE.
#' @param lambda.type specifies the value of lambda used for the predictions. \code{lambda.min} gives lambda with minimum cross-validated errors. \code{lambda.1se} gives the largest value of lambda such that the error is within 1 standard error of the minimum. Note that \code{lambda.1se} can only be chosen without restrictions of \code{max.coef}.
#' @param standardize logical, whether the predictors should be standardized or not. Default is TRUE.
#' @param nfolds the number of CV procedure folds.
#' @param foldid an optional vector of values between 1 and nfold identifying what fold each observation is in.
#' @param cvoffset logical, whether CV should be used to estimate the offsets. Default is FALSE.
#' @param cvoffsetnfolds the number of folds in the CV procedure that is performed to estimate the offsets. Default is 10. Only relevant if \code{cvoffset=TRUE}.
#' @param mcontrol controls how to deal with blockwise missing data. For details see below or \code{missing.control}
#' @param ... other arguments that can be passed to the function \code{cv.glmnet}.
#'
#' @return object of class \code{prioritylasso} with the following elements. If these elements are lists, they contain the results for each penalized block.
#' \describe{
#' \item{\code{lambda.ind}}{list with indices of lambda for \code{lambda.type}.}
#' \item{\code{lambda.type}}{type of lambda which is used for the predictions.}
#' \item{\code{lambda.min}}{list with values of lambda for \code{lambda.type}.}
#' \item{\code{min.cvm}}{list with the mean cross-validated errors for \code{lambda.type}.}
#' \item{\code{nzero}}{list with numbers of non-zero coefficients for \code{lambda.type}.}
#' \item{\code{glmnet.fit}}{list of fitted \code{glmnet} objects.}
#' \item{\code{name}}{a text string indicating type of measure.}
#' \item{\code{block1unpen}}{if \code{block1.penalization = FALSE}, the results of either the fitted \code{glm} or \code{coxph} object corresponding to \code{best.blocks}.}
#' \item{\code{coefficients}}{vector of estimated coefficients. If \code{block1.penalization = FALSE} and \code{family = gaussian} or \code{binomial}, the first entry contains an intercept.}
#' \item{\code{call}}{the function call.}
#' \item{\code{X}}{the original data used for the calculation}
#' \item{\code{missing.data}}{list with logical entries for every block which observation is missing (\code{TRUE} means missing)}
#' }
#'
#' @note The function description and the first example are based on the R package \code{ipflasso}. The second example is inspired by the example of \code{\link[glmnet]{cv.glmnet}} from the \code{glmnet} package.
#' @author Simon Klau, Roman Hornung, Alina Bauer \cr
#' Maintainer: Simon Klau (\email{simonklau@ibe.med.uni-muenchen.de})
#' @seealso \code{\link[prioritylasso]{pl_data}}, \code{\link[prioritylasso]{cvm_prioritylasso}}, \code{\link[ipflasso]{cvr.ipflasso}}, \code{\link[ipflasso]{cvr2.ipflasso}}
#' @references Klau, S., Jurinovic, V., Hornung, R., Herold, T., Boulesteix, A.-L. (2018). Priority-Lasso: a simple hierarchical approach to the prediction of clinical outcome using multi-omics data. BMC Bioinformatics 19, 322
#' @export
#' @import stats
#' @import glmnet
#' @import survival
#' @import utils
#' @importFrom checkmate assert_class
#' @examples
#' # gaussian
#'   prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian",
#'                 type.measure = "mse", blocks = list(bp1=1:75, bp2=76:200, bp3=201:500),
#'                 max.coef = c(Inf,8,5), block1.penalization = TRUE,
#'                 lambda.type = "lambda.min", standardize = TRUE, nfolds = 5, cvoffset = FALSE)
#'\dontrun{
#'   # cox
#'   # simulation of survival data:
#'   n <- 50;p <- 300
#'   nzc <- trunc(p/10)
#'   x <- matrix(rnorm(n*p), n, p)
#'   beta <- rnorm(nzc)
#'   fx <- x[, seq(nzc)]%*%beta/3
#'   hx <- exp(fx)
#'   # survival times:
#'   ty <- rexp(n,hx)
#'   # censoring indicator:
#'   tcens <- rbinom(n = n,prob = .3,size = 1)
#'   library(survival)
#'   y <- Surv(ty, 1-tcens)
#'   blocks <- list(bp1=1:20, bp2=21:200, bp3=201:300)
#'   # run prioritylasso:
#'   prioritylasso(x, y, family = "cox", type.measure = "deviance", blocks = blocks,
#'                 block1.penalization = TRUE, lambda.type = "lambda.min", standardize = TRUE,
#'                 nfolds = 5)
#'
#'   # binomial
#'   # using pl_data:
#'   prioritylasso(X = pl_data[,1:1028], Y = pl_data[,1029], family = "binomial", type.measure = "auc",
#'                 blocks = list(bp1=1:4, bp2=5:9, bp3=10:28, bp4=29:1028), standardize = FALSE)}
#'



prioritylasso <- function(X,
                          Y,
                          weights,
                          family = c("gaussian", "binomial", "cox"),
                          type.measure,
                          blocks,
                          max.coef = NULL,
                          block1.penalization = TRUE,
                          lambda.type = "lambda.min",
                          standardize = TRUE,
                          nfolds = 10,
                          foldid,
                          cvoffset = FALSE,
                          cvoffsetnfolds = 10,
                          mcontrol = missing.control(),
                          ...){
  
  
  if (packageVersion("glmnet") < "2.0.13") {
    stop("glmnet >= 2.0.13 needed for this function.", call. = FALSE)
  }
  
  if(is.null(max.coef)){
    max.coef <- rep(+Inf, length(blocks))
  } else {
    if(min(max.coef) < +Inf && lambda.type == "lambda.1se"){
      warning("lambda.1se can only be chosen without restrictions of max.coef and is set to lambda.min.")
      lambda.type = "lambda.min"
    }
    if (!setequal(length(blocks), length(max.coef))) {
      stop("The length of the entries of argument max.coef must equal the number of blocks.")
    }
  }
  
  
  if(sum(lapply(blocks, length) <= 1) != 0){
    stop("A block has to contain at least two predictors.")
  }
  
  if (anyDuplicated(as.numeric(unlist(blocks))) != 0 || !setequal(as.numeric(unlist(blocks)), 1:ncol(X))) {
    stop("Each predictor should be included in exactly one block.")
  }
  
  if (!is.element(lambda.type, c("lambda.min", "lambda.1se"))) {
    stop("lambda.type must be either lambda.min or lambda.1se.")
  }
  
  family <- match.arg(family)
  
  if (family == "gaussian") {
    if (type.measure != "mse")
      warning("type.measure is set to mse.")
    type.measure <- "mse"
  }
  if (family == "cox") {
    if (type.measure != "deviance")
      warning("type.measure is set to partial likelihood.")
    type.measure <- "deviance"
  }
  if (family == "binomial" & !is.element(type.measure, c("auc",
                                                         "class"))) {
    warning("type.measure is set to class.")
    type.measure <- "class"
  }
  
  if(type.measure == "auc") {
    if(cvoffset) {
      if(nrow(X)*((cvoffsetnfolds-1)/cvoffsetnfolds) - nrow(X)*((cvoffsetnfolds-1)/cvoffsetnfolds)*(nfolds-1)/nfolds < 10){
        stop(paste("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; Use nfolds < ", floor(nrow(X)*((cvoffsetnfolds-1)/cvoffsetnfolds)/10)+1, ".", sep=""))
      }
    }
    else {
      if(nrow(X) - nrow(X)*(nfolds-1)/nfolds < 10)
        stop(paste("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; Use nfolds < ", floor(nrow(X)/10)+1, ".", sep=""))
    }
  }
  
  
  if(missing(weights)){
    weights = rep(1, nrow(X))}
  else {
    if (length(weights) != nrow(X))
      stop(paste("number of elements in weights (", length(weights),
                 ") not equal to the number of rows of X (", nrow(X),
                 ")", sep = "")) }
  
  
  if(!missing(foldid)){
    if (length(foldid) != nrow(X))
      stop(paste("number of elements in foldid (", length(foldid),
                 ") not equal to the number of rows of X (", nrow(X),
                 ")", sep = ""))
    else {
      if(nfolds != max(foldid)){
        warning(paste("nfolds is set to", max(foldid)))
        nfolds = max(foldid)
      }
    }
  } else {
    foldid = sample(rep(seq(nfolds), length = nrow(X)))
  }
  
  # input check for handling missing data
  assert_class(mcontrol, "pl.missing.control")
  if (mcontrol$handle.missingdata == "none" && sum(is.na(X)) > 0) {
    stop("X contains missing data. Please use another value than 'none' for handle.missingdata.")
  }
  if (mcontrol$handle.missingdata != "none" && cvoffset) {
    stop("At the moment, a crossvalidated offset is only supported for complete data sets.")
  }
  
  # issue warnings concering missing data handling
  if (mcontrol$handle.missingdata == "ignore" ||
      mcontrol$handle.missingdata == "impute_offset") {
    foldid <- NULL
    warning(paste0("For handle.missingdata = ", mcontrol$handle.missingdata, ", the foldids of the observations are chosen individually for every block and not set globally. foldid is set to NULL"))
  }
  
  if (mcontrol$handle.missingdata == "impute_offset" &&
      mcontrol$impute.offset.cases == "complete.cases") {
    # calculate fraction of complete cases
    perc_complete_cases <- nrow(complete.cases(X)) / nrow(X)
    
    if (nrow(complete.cases(X)) == 0) {
      stop("The dataset contains no complete cases (over all blocks). Imputation of the offsets not possible.")
    }
    if (perc_complete_cases < mcontrol$perc.comp.cases.warning) {
      warning(paste0("The fraction of complete cases only is ",
                     round(perc_complete_cases, digits = 2)))
    }
    
    # check that every observation with missing data only has one missing block
    missing_index_overview <- matrix(FALSE, nrow = nrow(X),
                                     ncol = length(blocks))
    for (i in seq_along(blocks)) {
      missing_index_overview[, i] <- !complete.cases(X[, blocks[[i]]])
    }
    for (i in seq_along(nrow(missing_index_overview))) {
      if (sum(missing_index_overview[i, ]) > 1) {
        stop("For impute.offset.cases = 'complete.cases', every observation must only contain one missing block. ")
      }
    }
    
  }
  
  
  
  # generate a list which observations to use for which block
  # this is important for handle.missingdata = ignore
  observation_index <- lapply(blocks, function(block) {
    result <- which(complete.cases(X[, block]))
    if (length(result) == 0) {
      NULL
    } else {
      result
    }
  })
  missing_index <- lapply(blocks, function(block) {
    result <- which(!complete.cases(X[, block]))
    if (length(result) == 0) {
      NULL
    } else {
      result
    }
  })
  
  
  lambda.min <- list()
  lambda.ind <- list()
  min.cvm <- list()
  nzero <- list()
  glmnet.fit <- list()
  coeff <- list()
  lassoerg <- list()
  liste <- list(NULL)
  # list for every block, if TRUE the value is missing for this block
  missing.data <- list()
  start_block <- 1
  
  if (!block1.penalization) {
    if (length(blocks[[1]]) >= nrow(X)){
      stop("An unpenalized block 1 is only possible if the number of predictors in this block is smaller than the number of obervations.")
    }
    current_observations <- observation_index[[1]]
    current_missings <- missing_index[[1]]
    missing.data[[1]] <- seq(nrow(X)) %in% current_missings
    
    if (family != "cox") {
      block1erg <- glm(Y[current_observations] ~ X[current_observations,
                                                   blocks[[1]]],
                       family = family,
                       weights = weights[current_observations])
      predict_type <- "link"
      
    } else {
      block1erg <- coxph(Y[current_observations, ] ~ X[current_observations,
                                                       blocks[[1]]],
                         weights = weights[current_observations],
                         model = TRUE)
      predict_type <- "lp"
    }
    names(block1erg$coefficients) <- substr(names(block1erg$coefficients),
                                            start = 37,
                                            nchar(names(block1erg$coefficients)))
    
    if(cvoffset) {
      
      datablock1 <- data.frame(X[, blocks[[1]], drop = FALSE])
      datablock1$Y <- Y
      
      cvdiv <- makeCVdivision(n = nrow(X), K = cvoffsetnfolds, nrep = 1)[[1]]
      pred <- matrix(nrow = nrow(X), ncol = 1)
      for(count in seq(along = cvdiv)) {
        
        if (family != "cox") {
          block1ergtemp <- glm(Y ~ .,
                               data = datablock1[cvdiv[[count]] == 1, ],
                               weights = weights[cvdiv[[count]] == 1],
                               family = family)
          
          
          
        } else {
          block1ergtemp <- coxph(Y ~ ., data = datablock1[cvdiv[[count]] == 1, ],
                                 weights = weights[cvdiv[[count]] == 1])
        }
        
        names(block1ergtemp$coefficients) <- substr(names(block1ergtemp$coefficients),
                                                    start = 17,
                                                    nchar(names(block1ergtemp$coefficients)))
        pred[cvdiv[[count]]==0, ] <- as.matrix(predict(block1ergtemp,
                                                       newdata = datablock1[cvdiv[[count]] == 0, ]),
                                               type = predict_type)
        
      }
      
    } else {
      pred <- as.matrix(predict(block1erg, type = predict_type))
    }
    
    start_block <- 2
    #JH change here or add check to forbid missing data in first block until
    # further clarification
    # calculate the new offsets
    new_offsets <- calculate_offsets(current_missings = current_missings,
                                     current_observations = current_observations,
                                     mcontrol = mcontrol,
                                     current_block = 1,
                                     pred = pred,
                                     liste = liste,
                                     X = X,
                                     blocks = blocks)
    liste[[2]] <- new_offsets
    lassoerg <- list(block1erg)
    coeff[[1]] <- block1erg$coefficients
  } else {
    block1erg <- NULL
  }
  
  
  for (i in start_block:length(blocks)) {
    
    actual_block <- blocks[[i]]
    current_observations <- observation_index[[i]]
    current_missings <- missing_index[[i]]
    missing.data[[i]] <- seq(nrow(X)) %in% current_missings
    
    # if i == 1 for the first block, then offset is NULL, because first list
    # element is NULL
    
    # when foldid <- NULL for handle.missingdata  != "none",
    # foldid[current_observations] still evaluates to NULL
    lassoerg[[i]] <- cv.glmnet(X[current_observations, actual_block],
                               Y[current_observations],
                               weights[current_observations],
                               offset = liste[[i]][current_observations],
                               family = family,
                               type.measure = type.measure,
                               nfolds = nfolds,
                               foldid = foldid[current_observations],
                               standardize = standardize,
                               ...)
    
    # determine the lambdas with the best results (and which are used for the
    # predictions if no crossvalidation of the offsets are used)
    if (lambda.type == "lambda.1se") {
      lambda_to_use <- "lambda.1se"
      
      lambda.ind[i] <- which(lassoerg[[i]]$lambda == lassoerg[[i]][lambda.type])
      lambda.min[i] <- lassoerg[[i]][lambda.type]
      
    } else {
      # if lambda.type == "lambda.min" calculate the correct lambda
      # with respect to max.coeff (if not specified, max.coef ist +Inf)
      which_lambda <- which(as.numeric(lassoerg[[i]]$nzero) <= max.coef[i])
      
      if (type.measure != "auc"){
        lcvmi <- lassoerg[[i]]$cvm
      } else {
        lcvmi <- -lassoerg[[i]]$cvm
      }
      
      lambda_to_use <- lassoerg[[i]]$lambda[which_lambda[which.min(lcvmi[which_lambda])[1]]]
      
      lambda.min[i] <- lambda_to_use
      lambda.ind[i] <- which(lassoerg[[i]]$lambda == lambda.min[i])
    }
    
    
    
    if(cvoffset) {
      cvdiv <- makeCVdivision(n = nrow(X), K = cvoffsetnfolds, nrep = 1)[[1]]
      pred <- matrix(nrow = nrow(X), ncol = 1)
      for(count in seq(along=cvdiv)) {
        lassoergtemp <- cv.glmnet(X[cvdiv[[count]] == 1,actual_block, drop = FALSE],
                                  Y[cvdiv[[count]] == 1],
                                  weights[cvdiv[[count]] == 1],
                                  offset = liste[[i]][cvdiv[[count]] == 1, drop = FALSE],
                                  family = family,
                                  type.measure = type.measure,
                                  nfolds = nfolds,
                                  foldid = foldid[cvdiv[[count]] == 1],
                                  standardize = standardize,
                                  ...)
        
        # determine the correct lambda to use for the predictions
        if (lambda.type == "lambda.1se") {
          lambda_to_use <- "lambda.1se"
        } else {
          # if lambda.type == "lambda.min" calculate the correct lambda
          # with respect to max.coeff (if not specified, max.coeff is +Inf)
          which_lambdatemp <- which(as.numeric(lassoergtemp$nzero) <=
                                      max.coef[i])
          
          if (type.measure != "auc"){
            lcvmitemp <-  lassoergtemp$cvm
          } else {
            lcvmitemp <-  -lassoergtemp$cvm
          }
          
          lambda_to_use <- lassoergtemp$lambda[which_lambdatemp[which.min(lcvmitemp[which_lambdatemp])[1]]]
          
        }
        
        
        pred[cvdiv[[count]] == 0, ] <- predict(lassoergtemp,
                                               newx = X[cvdiv[[count]] == 0, actual_block],
                                               newoffset = liste[[i]][cvdiv[[count]] == 1, drop = FALSE],
                                               s = lambda_to_use,
                                               type = "link")
        
      }
    }
    else {
      pred <- predict(lassoerg[[i]],
                      newx = X[current_observations, actual_block],
                      newoffset = liste[[i]][current_observations],
                      s = lambda_to_use,
                      type = "link")
    }
    
    # store the results for the current block
    # calculate the new offsets
    new_offsets <- calculate_offsets(current_missings = current_missings,
                                     current_observations = current_observations,
                                     mcontrol = mcontrol,
                                     current_block = i,
                                     pred = pred,
                                     liste = liste,
                                     X = X,
                                     blocks = blocks)
    liste[[i+1]] <- new_offsets
    
    min.cvm[i] <- lassoerg[[i]]$cvm[lambda.ind[[i]]]
    nzero[i] <- lassoerg[[i]]$nzero[lambda.ind[[i]]]
    glmnet.fit[[i]] <- lassoerg[[i]]$glmnet.fit
    coeff[[i]] <- glmnet.fit[[i]]$beta[,lambda.ind[[i]]]
  }
  
  
  name <- lassoerg[[i]]$name
  
  
  
  
  finallist <- list(lambda.ind = lambda.ind,
                    lambda.type = lambda.type,
                    lambda.min = lambda.min,
                    min.cvm = min.cvm,
                    nzero = nzero,
                    glmnet.fit = glmnet.fit,
                    name = name,
                    block1unpen = block1erg,
                    coefficients = unlist(coeff),
                    call = match.call(),
                    X = X,
                    missing.data = missing.data)
  
  class(finallist) <- c("prioritylasso", class(finallist))
  return(finallist)
  
}
