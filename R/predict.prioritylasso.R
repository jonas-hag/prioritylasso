#' Predictions from prioritylasso
#'
#' Makes predictions for a \code{prioritylasso} object. It can be chosen between linear predictors or fitted values.
#'
#' \code{handle.missingtestdata} specifies how to deal with missing data.
#' The default \code{none} cannot handle missing data, \code{set.zero} leaves out
#' the missing data for the calculation of the prediction (the missing value is set to zero).
#' \code{impute.block} uses an imputation model to impute the offset of a missing block. This only works if the prioritylasso object was fitted with \code{handle.missingdata = "impute.offset"}.
#' If \code{impute.offset.cases = "complete.cases"} was used, then every observation can have only one missing block. If \code{impute.offset.cases = "available.cases"} was used, the missingness pattern in the test data has to be the same as in the train data.
#'
#' @param object An object of class \code{prioritylasso}.
#' @param newdata (nnew \code{x} p) matrix or data frame with new values.
#' @param type Specifies the type of predictions. \code{link} gives the linear predictors for all types of response and \code{response} gives the fitted values.
#' @param handle.missingtestdata Specifies how to deal with missing data in the test data; possibilities are \code{none}, \code{set.zero} and \code{impute.block}
#' @param include.allintercepts should the intercepts from all blocks included in the prediction? If \code{FALSE}, only the intercept from the first block is included (default in the past).
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Predictions that depend on \code{type}.
#'
#' @author Simon Klau
#' @export
#' @importFrom checkmate assert_logical
#' @import glmnet
#' @seealso \code{\link[prioritylasso]{pl_data}}, \code{\link[prioritylasso]{prioritylasso}}
#' @examples
#' pl_bin <- prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rbinom(50,1,0.5),
#'                        family = "binomial", type.measure = "auc",
#'                        blocks = list(block1=1:13,block2=14:200, block3=201:500),
#'                        block1.penalization = TRUE, lambda.type = "lambda.min",
#'                        standardize = FALSE, nfolds = 5)
#'
#' newdata_bin <- matrix(rnorm(30*500),30,500)
#'
#' predict(object = pl_bin, newdata = newdata_bin, type = "response")




predict.prioritylasso <- function(object,
                                  newdata = NULL,
                                  type = c("link", "response"),
                                  handle.missingtestdata = c("none", "set.zero",
                                                             "impute.block"),
                                  include.allintercepts = FALSE,
                                  ...){
  
  # input check
  type <- match.arg(type)
  handle.missingtestdata <- match.arg(handle.missingtestdata)
  assert_logical(include.allintercepts)
  
  if (is.null(newdata)) {
    newdata <- as.matrix(object$X)
  } else {
    newdata <- data.matrix(newdata)
  }
  
  if (ncol(newdata) != ncol(object$X)) {
    stop("The newdata does not have the same number of covariates as the traindata in the prioritylasso object.")
  }
  
  
  if (handle.missingtestdata == "none" && sum(is.na(newdata)) > 0) {
    stop("X contains missing data. Please use another value than 'none' for handle.missingtestdata.")
  }
  
  # check that there are no single missing values
  lapply(object$blocks, function(block) {
    lapply(seq_len(nrow(newdata)), function(i) {
      if (sum(is.na(newdata[i, block])) != 0 &&
          sum(is.na(newdata[i, block])) != length(block)) {
        stop(paste0("Observation ", i, " contains a single missing value. This is not supported."))
      }
    })
  })
  
  # check that there are no observations with no values at all
  lapply(seq_len(nrow(newdata)), function(i) {
    if (sum(is.na(newdata[i, ])) == ncol(newdata)) {
      stop(paste0("Observation ", i, " only consists of missing values."))
    }
  })
  
  # check that if there are only complete cases, handle.missingtestdata = none
  if (handle.missingtestdata != "none" &&
      sum(complete.cases(newdata)) == nrow(newdata)) {
    stop("The data consists only of complete cases. Please use handle.missingtestdata = 'none'.")
  }
  
  if (handle.missingtestdata == "impute.block") {
    # check that the prioritylasso object contains the imputation models
    if (object$call$mcontrol$handle.missingdata != "impute.offset" ||
        is.null(object$imputation.models)) {
      stop("The needed imputation models are not provided by the prioritylasso object. Refit the training data with handle.missingdata = 'impute.offset")
    }
    
    
    # if the predicted value of a missing block should be imputed, check that
    #  - there is only one missing block per observation (if the original model
    # was fitted with impute.offset.cases = complete.cases)
    #  - there is only the same missingness pattern as in the training data
    # (if the original model was fitted with impute.offset.cases = available.cases)
    missing_index_overview <- matrix(FALSE, nrow = nrow(newdata),
                                     ncol = length(object$blocks))
    for (i in seq_along(object$blocks)) {
      missing_index_overview[, i] <- !complete.cases(newdata[, object$blocks[[i]]])
    }
    
    # for complete.cases (is standard in missing.control, so it can be NULL in
    # the call)
    if (is.null(object$call$mcontrol$impute.offset.cases) ||
        object$call$mcontrol$impute.offset.cases == "complete.cases") {
      for (i in seq_len(nrow(missing_index_overview))) {
        if (sum(missing_index_overview[i, ]) > 1) {
          stop("For handle.missingtestdata = 'impute.block', every observation must only contain one missing block.")
        }
      }
    } else {
      # for available.cases
      # for every observation, check if a block is missing. If this is the case,
      # check that the required blocks for the imputation model have no missing
      # data
      # for every observation, check if a block is missing. If this is the case,
      # check that observation fits to a missingness pattern observed in the
      # training data
      missing_pattern_overview <- matrix(NA, nrow = nrow(newdata),
                                       ncol = length(object$blocks))
      for (i in seq_len(nrow(missing_index_overview))) {
        for (j in seq_len(ncol(missing_index_overview))) {
          # TRUE means that the block is missing
          if (missing_index_overview[i, j]) {
            # if a block j is missing, check that the blocks required to impute
            # this block j have no missing data
            # if no missingness pattern from the training data fits the
            # missingness pattern of the current observation, the sum of the
            # comparison is 0
            if (sum(compare_boolean(object$missingness.pattern[[j]],
                                    missing_index_overview[i, ])) == 0) {
              stop(paste0("Observation ", i, " has the missingness pattern ",
                          paste0(missing_index_overview[i, ], collapse = " "),
                          ". This pattern was not present in the training data, therefore no value can be predicted."))
            } else {
              # determine the correct missingness pattern (-> and attached
              # imputation model) for every observation/block combination
              missing_pattern_overview[i, j] <-
                which(compare_boolean(object$missingness.pattern[[j]],
                                      missing_index_overview[i, ]))
            }
          }
        }
      }
      
      
    }
  }
  
  # if blocks should be imputed -> check for every observation if/which blocks
  # needs imputation
  # -> calculate this imputation and store it in add_vector (with 0 for non
  # imputed observations); this vector is added to the prediction
  #
  # set the NAs for the imputed block in this observation to 0
  
  # for every observation, make a matrix with 1 if this intercept should be
  # included or 0 if not (needed for imputed observations that the intercept
  # is not added additionally)
  
  # generate vector with intercepts of the different blocks
  # cox model does not have intercepts
  if (object$call$family == "cox") {
    intercepts <- rep(0, times = length(object$blocks))
  } else {
    # unpenalized first block
    if (!is.null(object$block1unpen)) {
      intercepts <- object$coefficients[1]
      for (i in 2:length(object$glmnet.fit)) {
        intercepts[i] <- object$glmnet.fit[[i]]$a0[object$lambda.ind[[i]]]
      }
    } else {
      # penalized first block
      intercepts <- c()
      for (i in 1:length(object$glmnet.fit)) {
        intercepts[i] <- object$glmnet.fit[[i]]$a0[object$lambda.ind[[i]]]
      }
    }
  }
  
  # generate matrix with entry for every observation (row) and every intercept
  # (column) if it should be included
  intercept_model_matrix <- matrix(0, nrow = nrow(newdata),
                                   ncol = length(object$blocks))
  
  # determine if only the first intercept or all intercepts should be included
  if (include.allintercepts) {
    intercept_model_matrix[, ] <- 1
  } else {
    intercept_model_matrix[, 1] <- 1
  }
  
  imputed_values <- rep(0, nrow(newdata))
  if (handle.missingtestdata == "impute.block") {
    # determine for which observation a block has to be imputed
    # determine which block (or which not -> NA) has to be imputed for every
    # observation
    # there can be several blocks per observation that get imputed
    impute_which_block <- apply(missing_index_overview, 1, function(x) {
      res <- which(x)
      if (length(res) == 0) {
        NA
      } else {
        res
      }
    })
    # determine for which observations a block has to be imputed
    index_observation <- which(!is.na(impute_which_block))
    impute_which_block <- impute_which_block[index_observation]
    missing_pattern_overview <- missing_pattern_overview[index_observation, ]
    
    # for these observations (i) and the block (j), the corresponding intercept
    # doesn't need to be included
    for (i in 1:length(index_observation)) {
      for (j in seq_along(impute_which_block[[i]])) {
        intercept_model_matrix[index_observation[i],
                               impute_which_block[[i]][j]] <- 0
      }
    }
    
    # perform the imputations for the missing data
    if (is.null(object$call$mcontrol$lambda.imputation)) {
      use_lambda <- missing.control()[["lambda.imputation"]]
    } else {
      use_lambda <- object$call$mcontrol$lambda.imputation
    }
    for (i in 1:length(index_observation)) {
      # for every observation, check the different blocks
      ################ explain the indexing
      # index_observation and impute_which_block contain already only these
      # observations where something has to be imputed
      # therefore, i gives the index starting from 1 with which observation is
      # dealt with, and j gives the missing block (absolute, no index)
      # impute_which_block[[i]] gives the vector which blocks have to be
      # imputed for the current observation
      # impute_which_block[[i]][j] gives one block of this observation
      # missing_pattern_overview[i, j] gives the information which of the
      # possible several imputation models per block has to be used
      # -> with this, you can chose the imputation model for the corresponding 
      # block
      # -> also, you can chose for this block which other blocks are used in
      # this imputation model
      # with this, you can select the corresponding columns for these blocks
      # because they are returned as a list, unlist them
      for (j in impute_which_block[[i]]) {
        # check if a model exists or if it is just a constant
        if (class(object$imputation.models[[j]][[missing_pattern_overview[i, j]]]) == "cv.glmnet") {
          imputed_values[index_observation[i]] <-
            predict(object$imputation.models[[j]][[missing_pattern_overview[i, j]]],
                    newx = newdata[index_observation[i],
                                   unlist(object$blocks[object$blocks.used.for.imputation[[j]][[missing_pattern_overview[i, j]]]]),
                                   drop = FALSE],
                    s = use_lambda) +
            imputed_values[index_observation[i]]
          
        } else {
          if (class(object$imputation.models[[j]][[missing_pattern_overview[i, j]]]) == "numeric") {
            imputed_values[index_observation[i]] <-
              object$imputation.models[[j]][[missing_pattern_overview[i, j]]] +
              imputed_values[index_observation[i]]
          } else {
            imputed_values[index_observation[i]] <- NA 
          }
        }
      }
    }
  }
  
  # calculate the (several) intercept values for every observation
  add_intercept <- intercept_model_matrix %*% intercepts
  
  # if missing data should be ignored, set NAs to 0
  # the same if blocks should get imputed
  if (handle.missingtestdata == "set.zero" ||
      handle.missingtestdata == "impute.block") {
    newdata <- replace(newdata, which(is.na(newdata)), 0)
  }
  
  # determine the coefficients
  # family = cox
  if (object$call$family == "cox") {
    
    coeff <- object$coefficients
    
  } else {
    # family = gaussian/binomial
    # first block is penalized
    if (is.null(object$block1unpen)) {
      
      coeff <- object$coefficients
      
      # family = gaussian/binomial
      # first block is not penalized
    } else {
      
      coeff <- object$coefficients[-1]
    }
  }
  
  pred <- newdata %*% coeff + add_intercept + imputed_values
  
  
  if(type == "response"){
    if(object$call$family == "binomial"){
      pred <- exp(pred) / (1 + exp(pred)) # fitted probabilities
    }
    if(object$call$family == "cox"){
      pred <- exp(pred) # fitted relative risk (risk score (exp(lp)))
    }
  }
  
  # rescaling if gaussian outcome was normalised in training data
  if (object$call$family == "gaussian" &&
      !is.null(object$y.scale.param)) {
    pred <- pred * object$y.scale.param$sd + object$y.scale.param$mean
  }
  return(pred)
}
