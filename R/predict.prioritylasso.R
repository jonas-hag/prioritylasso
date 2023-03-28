#' Predictions from prioritylasso
#'
#' Makes predictions for a \code{prioritylasso} object. It can be chosen between linear predictors or fitted values.
#'
#' \code{handle.missingtestdata} specifies how to deal with missing data.
#' The default \code{none} cannot handle missing data, \code{omit.prediction} does not make a prediction for observations with missing values and return \code{NA}. \code{set.zero} ignores
#' the missing data for the calculation of the prediction (the missing value is set to zero).
#' \code{impute.block} uses an imputation model to impute the offset of a missing block. This only works if the prioritylasso object was fitted with \code{handle.missingdata = "impute.offset"}.
#' If \code{impute.offset.cases = "complete.cases"} was used, then every observation can have only one missing block. For observations with more than one missing block, \code{NA} is returned.
#' If \code{impute.offset.cases = "available.cases"} was used, the missingness pattern in the test data has to be the same as in the train data. For observations with an unknown missingness pattern, \code{NA} is returned.
#'
#' @param object An object of class \code{prioritylasso}.
#' @param newdata (nnew \code{x} p) matrix or data frame with new values.
#' @param type Specifies the type of predictions. \code{link} gives the linear predictors for all types of response and \code{response} gives the fitted values.
#' @param handle.missingtestdata Specifies how to deal with missing data in the test data; possibilities are \code{none}, \code{omit.prediction}, \code{set.zero} and \code{impute.block}
#' @param include.allintercepts should the intercepts from all blocks included in the prediction? If \code{FALSE}, only the intercept from the first block is included (default in the past).
#' @param use.blocks determines which blocks are used for the prediction, the default is all. Otherwise one can specify the number of blocks which are used in a vector
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Predictions that depend on \code{type}.
#'
#' @author Simon Klau
#' @export
#' @importFrom checkmate assert_logical assert_numeric assert check_matrix check_data_frame
#' @import glmnet
#' @seealso \code{\link[prioritylasso]{pl_data}}, \code{\link[prioritylasso]{prioritylasso}}
#' @examples
#' pl_bin <- prioritylasso(X = matrix(rnorm(50*190),50,190), Y = rbinom(50,1,0.5),
#'                        family = "binomial", type.measure = "auc",
#'                        blocks = list(block1=1:13,block2=14:80, block3=81:190),
#'                        block1.penalization = TRUE, lambda.type = "lambda.min",
#'                        standardize = FALSE, nfolds = 3)
#'
#' newdata_bin <- matrix(rnorm(10*190),10,190)
#'
#' predict(object = pl_bin, newdata = newdata_bin, type = "response")




predict.prioritylasso <- function(object,
                                  newdata = NULL,
                                  type = c("link", "response"),
                                  handle.missingtestdata = c("none",
                                                             "omit.prediction",
                                                             "set.zero",
                                                             "impute.block"),
                                  include.allintercepts = FALSE,
                                  use.blocks = "all",
                                  ...){
  
  ##############################################################################
  ############################ input check
  ##############################################################################
  type <- match.arg(type)
  handle.missingtestdata <- match.arg(handle.missingtestdata)
  assert_logical(include.allintercepts)
  if (use.blocks[1] != "all" || length(use.blocks) != 1) {
    # if not all blocks are selected via "all", one can only select at max. all
    # blocks with their order from 1 to last block
    assert_numeric(use.blocks, min.len = 1, max.len = length(object$blocks),
                   unique = TRUE, lower = 1, upper = length(object$blocks))
  }
  
  if (is.null(newdata)) {
    if (!inherits(object$X, "matrix") || all(is.na(object$X))) {
      stop("No data provided by either the prioritylasso object or newdata.")
    } else {
      newdata <- as.matrix(object$X)
    }
  } else {
    assert(check_matrix(newdata), check_data_frame(newdata))
    newdata <- data.matrix(newdata)
  }
  
  if (ncol(newdata) != object$dim.x[2]) {
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
  
  # store one parameter from the model in a local variable for easier access
  if (is.null(object$mcontrol$impute.offset.cases)) {
    # this is the default from missing.control
    param_offset_cases <- "complete.cases"
  } else {
    param_offset_cases <- object$mcontrol$impute.offset.cases
  }
  
  ##############################################################################
  ############################ omit.prediction
  ############################ exclude observations
  ##############################################################################
  
  # if handle.missingtestdata = omit.prediction, search for observations with
  # missing values -> these are excluded from the prediction
  index_exclude_observations <- NULL
  missing_index_overview <- matrix(FALSE, nrow = nrow(newdata),
                                   ncol = length(object$blocks))
  if (handle.missingtestdata == "omit.prediction") {
    for (i in seq_len(nrow(missing_index_overview))) {
      if (sum(missing_index_overview[i, ]) > 0) {
        index_exclude_observations <- c(index_exclude_observations, i)
      }
    }
    
    # exclude observations with missing values
    index_include_observations <- setdiff(seq_len(nrow(newdata)),
                                          index_exclude_observations)
    newdata <- newdata[index_include_observations, ]
    missing_index_overview <-
      missing_index_overview[index_include_observations, ]
  }
  
  ##############################################################################
  ############################ impute.block
  ############################ exclude observations
  ##############################################################################
  
  # generate the index which blocks are missing for which observation
  # (outside of conditions only for impute.block, as the index is also needed
  # to exclude intercepts from the prediction for handle.missingtestdata =
  # set.zero when the model was trained with handle.missingdata = ignore)
  for (i in seq_along(object$blocks)) {
    missing_index_overview[, i] <- !complete.cases(newdata[, object$blocks[[i]]])
  }
  
  if (handle.missingtestdata == "impute.block") {
    # check that the prioritylasso object contains the imputation models
    if (object$mcontrol$handle.missingdata != "impute.offset" ||
        is.null(object$imputation.models)) {
      stop("The needed imputation models are not provided by the prioritylasso object. Refit the training data with handle.missingdata = 'impute.offset")
    }
    
    # if handle.missingtestdata = impute.block,
    # all observations with unfitting missingness patterns are excluded and NA
    # returned for this observations
    
    # if the predicted value of a missing block should be imputed, check that
    #  - there is only one missing block per observation (if the original model
    # was fitted with impute.offset.cases = complete.cases)
    #  - there is only the same missingness pattern as in the training data
    # (if the original model was fitted with impute.offset.cases = available.cases)
    
    
    # for complete.cases
    if (param_offset_cases == "complete.cases") {
      for (i in seq_len(nrow(missing_index_overview))) {
        if (sum(missing_index_overview[i, ]) > 1) {
          index_exclude_observations <- c(index_exclude_observations, i)
          # stop("For handle.missingtestdata = 'impute.block', every observation must only contain one missing block.")
        }
      }
      
      index_include_observations <- setdiff(seq_len(nrow(newdata)),
                                            index_exclude_observations)
      
    } else {
      
      ##########################################################################
      ############################ available.cases
      ############################ determine missingness pattern
      ##########################################################################
      
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
              index_exclude_observations <- c(index_exclude_observations, i)
              
              warning(paste0("Observation ", i, " has the missingness pattern ",
                          paste0(missing_index_overview[i, ], collapse = " "),
                          ". This pattern was not present in the training data, therefore no value can be predicted for this observation."))
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
      # due to checking every block separately, the index could contain an
      # observation several times
      index_exclude_observations <- unique(index_exclude_observations)
      
      # delete the rows in missing_pattern_overview for excluded observations
      index_include_observations <- setdiff(seq_len(nrow(newdata)),
                                            index_exclude_observations)
      missing_pattern_overview <-
        missing_pattern_overview[index_include_observations, ]
      
    }
    # exclude observations with missing values
    newdata <- newdata[index_include_observations, ]
    missing_index_overview <-
      missing_index_overview[index_include_observations, , drop = FALSE]
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
  
  ##############################################################################
  ############################ get intercepts
  ##############################################################################
  
  # generate vector with intercepts of the different blocks
  # cox model does not have intercepts
  if (object$family == "cox") {
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
  
  ##############################################################################
  ############################ determine which intercepts to use
  ##############################################################################
  
  # generate matrix with entry for every observation (row) and every intercept
  # (column) if it should be included
  intercept_model_matrix <- matrix(0, nrow = nrow(newdata),
                                   ncol = length(object$blocks))
  
  # determine if only the first intercept or all intercepts should be included
  # also take into account that only these intercepts are used that belong to
  # to the blocks which should be used for prediction
  if (include.allintercepts) {
    # only use the intercepts from the blocks used for prediction
    if (use.blocks[1] == "all") {
      block_index <- 1:length(object$blocks)
    } else {
      block_index <- use.blocks
    }
    
    intercept_model_matrix[, block_index] <- 1
    
  } else {
    if (use.blocks[1] == "all" || 1 %in% use.blocks) {
      intercept_model_matrix[, 1] <- 1
    }
  }
  
  # if handle.missingtestdata = set.zero, also exclude the intercept of the
  # blocks where the data is set to zero
  if (handle.missingtestdata == "set.zero") {
    # check for every observation
    for (i in seq_len(nrow(intercept_model_matrix))) {
      # for every block
      for (j in seq_len(ncol(intercept_model_matrix))) {
        # if the observation is missing this block, also exclude the intercept
        if (missing_index_overview[i, j]) {
          intercept_model_matrix[i, j] <- 0
        }
      }
    }
  }
  
  ##############################################################################
  ############################ do the imputation if necessary
  ##############################################################################
  
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
    
    # if the missing data in testdata is completely removed because the pattern
    # is not in the training data, then no observations have to be imputed
    # (and the length of the index is 0)
    if (length(index_observation) != 0) {
      impute_which_block <- impute_which_block[index_observation]
      if (param_offset_cases == "available.cases") {
        missing_pattern_overview <- missing_pattern_overview[index_observation, , drop = FALSE]
      }
      
      # for these observations (i) and the block (j), the corresponding intercept
      # doesn't need to be included
      for (i in 1:length(index_observation)) {
        for (j in seq_along(impute_which_block[[i]])) {
          intercept_model_matrix[index_observation[i],
                                 impute_which_block[[i]][j]] <- 0
        }
      }
      
      # perform the imputations for the missing data
      if (is.null(object$mcontrol$lambda.imputation)) {
        use_lambda <- missing.control()[["lambda.imputation"]]
      } else {
        use_lambda <- object$mcontrol$lambda.imputation
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
        
        ########################################################################
        ############################ imputation for available.cases
        ########################################################################
        if (param_offset_cases == "available.cases") {
          for (j in impute_which_block[[i]]) {
            # check if the current block should be used for prediction
            if (j %in% use.blocks || use.blocks[1] == "all") {
              # check if a model exists or if it is just a constant
              if (inherits(object$imputation.models[[j]][[missing_pattern_overview[i, j]]], "cv.glmnet")) {
                imputed_values[index_observation[i]] <-
                  predict(object$imputation.models[[j]][[missing_pattern_overview[i, j]]],
                          newx = newdata[index_observation[i],
                                         unlist(object$blocks[object$blocks.used.for.imputation[[j]][[missing_pattern_overview[i, j]]]]),
                                         drop = FALSE],
                          s = use_lambda) +
                  imputed_values[index_observation[i]]
                
              } else {
                if (inherits(object$imputation.models[[j]][[missing_pattern_overview[i, j]]], "numeric")) {
                  imputed_values[index_observation[i]] <-
                    object$imputation.models[[j]][[missing_pattern_overview[i, j]]] +
                    imputed_values[index_observation[i]]
                } else {
                  imputed_values[index_observation[i]] <- NA 
                }
              }
            }
            
          }
        } else {
          ######################################################################
          ############################ imputation for complete.cases
          ######################################################################
          imputed_block <- impute_which_block[[i]]
          if (imputed_block %in% use.blocks || use.blocks[1] == "all") {
            # check if a model exists or if it is just a constant
            if (inherits(object$imputation.models[[imputed_block]], "cv.glmnet")) {
              imputed_values[index_observation[i]] <-
                predict(object$imputation.models[[imputed_block]],
                        newx = newdata[index_observation[i],
                                       unlist(object$blocks[object$blocks.used.for.imputation[[imputed_block]]]),
                                       drop = FALSE],
                        s = use_lambda) +
                imputed_values[index_observation[i]]
              
            } else {
              if (inherits(object$imputation.models[[imputed_block]], "numeric")) {
                imputed_values[index_observation[i]] <-
                  object$imputation.models[[imputed_block]] +
                  imputed_values[index_observation[i]]
              } else {
                imputed_values[index_observation[i]] <- NA 
              }
            }
          }
          
        }
      }
    }
  }
  
  # calculate the (several) intercept values for every observation
  add_intercept <- intercept_model_matrix %*% intercepts
  
  # if missing data should be ignored, set NAs to 0
  # the same if blocks should get imputed (because the blocks with NA were
  # already imputed in the step above)
  if (handle.missingtestdata == "set.zero" ||
      handle.missingtestdata == "impute.block") {
    newdata <- replace(newdata, which(is.na(newdata)), 0)
  }
  
  ##############################################################################
  ############################ determine the coefficients
  ##############################################################################
  # family = cox
  if (object$family == "cox") {
    
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
  
  # only take the coefficients and covariates of the blocks that should be used
  # for prediction
  if (use.blocks[1] != "all") {
    index_covariates <- unlist(lapply(use.blocks, function(i) {
      object$blocks[[i]]
    }))
    index_covariates <- sort(index_covariates)
    coeff <- coeff[index_covariates]
    newdata <- newdata[, index_covariates]
  }
  
  ##############################################################################
  ############################ calculate the prediction & finetuning
  ##############################################################################
  
  pred <- newdata %*% coeff + add_intercept + imputed_values
  
  
  if(type == "response"){
    if(object$family == "binomial"){
      pred <- exp(pred) / (1 + exp(pred)) # fitted probabilities
    }
    if(object$family == "cox"){
      pred <- exp(pred) # fitted relative risk (risk score (exp(lp)))
    }
  }
  
  # rescaling if gaussian outcome was normalised in training data
  if (object$family == "gaussian" &&
      !is.null(object$y.scale.param)) {
    pred <- pred * object$y.scale.param$sd + object$y.scale.param$mean
  }
  
  # if for some obsevations no prediction could be made (due to missing values),
  # return NA for these observations
  if (length(index_exclude_observations) > 0) {
    # add the position index of the predicted observations
    pred <- cbind(pred, index_include_observations)
    
    # add NAs for the observations, for which a prediction is not possible
    na_matrix <- matrix(NA, nrow = length(index_exclude_observations),
                        ncol = 1)
    na_matrix <- cbind(na_matrix, index_exclude_observations)
    pred <- rbind(pred, na_matrix)
    
    # bring the observations into the correct order
    pred <- pred[order(pred[, 2]), ]
    pred <- pred[, 1, drop = FALSE]
    
  }
  
  return(pred)
}
