#' Calculates the offsets for the current block
#'
#' @param current_missings index vector (indices) of current missing observations
#' @param current_observations index vector (indices) of current used observations
#' @param mcontrol control for missing data handling
#' @param current_block index of current block
#' @param pred predictions of current block
#' @param liste list with offsets
#' @param X original data
#' @param blocks information which variable belongs to which block
#' @param current_intercept the intercept estimated for the current block
#'
#' @return list with offsets, used imputation model and the blocks used for the imputation model (if applicable)
calculate_offsets <- function(current_missings,
                              current_observations,
                              mcontrol,
                              current_block,
                              pred,
                              liste,
                              X,
                              blocks,
                              current_intercept) {
  
  # store the results for the current block
  # calculate the new offsets
  imputation_model <- NULL
  blocks_used_for_imputation <- NULL
  # for handle.missingdata = impute.offset, an imputation model has to be
  # calculated always, even if there are no missing values
  if (is.null(current_missings) &&
      mcontrol$handle.missingdata != "impute.offset") {
    new_offsets <- as.matrix(pred)
  } else {
    # the calculated offsets for this block where there are observations
    calculated_offsets <- as.matrix(pred)
    calculated_offsets <- cbind(calculated_offsets, current_observations)
    missing_offsets <- NULL
    # if chosen, ignore the missing block
    if (mcontrol$handle.missingdata == "ignore") {
      # for the missing values, take the offset from the previous block
      # if the missings are in the first block, use 0 as offset
      #JH -> needs further clarifaction
      if (current_block == 1) {
        if (mcontrol$offset.firstblock == "zero") {
          missing_offsets <- cbind(rep(0, length(current_missings)),
                                   current_missings)
        } else {
          missing_offsets <- cbind(rep(current_intercept,
                                       length(current_missings)),
                                   current_missings)
        }
      } else {
        missing_offsets <- cbind(liste[[current_block]][current_missings],
                                 current_missings)
      }
    }
    # if chosen impute the missing offsets
    if (mcontrol$handle.missingdata == "impute.offset") {
      # mcontrol$impute.offset.cases can either be complete.cases or
      # available.cases, has to be checked by input check beforehand
      if (mcontrol$impute.offset.cases == "complete.cases") {
        # get the x values for the imputation
        x_values <- cbind(X, 1:nrow(X))
        # only take the complete cases over all blocks
        x_values <- x_values[complete.cases(x_values), ]
        # for this imputation model, exclude the information about the current
        # block (because this would completeley determine the offset)
        x_values <- x_values[, -blocks[[current_block]]]
        blocks_used_for_imputation <- setdiff(seq_along(blocks), current_block)
      }
      
      if (mcontrol$impute.offset.cases == "available.cases") {
        # get the x values for the imputation
        x_values <- cbind(X, row_index = 1:nrow(X))
        # only take the observations which have values for the current block
        # and other blocks
        
        # get all the x values that have information for the current block
        row_index <- complete.cases(x_values[, blocks[[current_block]]])
        x_values <- x_values[row_index, ]
        
        # get further information which blocks can be used for the imputation
        # model
        # -> get all the blocks for which the current missing values have
        # information
        na_matrix_current_missings <- is.na(X[current_missings, ])
        blocks_without_missing <- c()
        blocks_to_check <- setdiff(seq_along(blocks), current_block)
        # check for every other than the current block if all the current
        # missing observations have information in the corresponding block
        # if yes, it can be used for the imputation model
        for (i in blocks_to_check) {
          if (sum(na_matrix_current_missings[, blocks[[i]]]) == 0) {
            blocks_without_missing <- c(blocks_without_missing, i)
          }
        }
        
        if (length(blocks_without_missing) == 0) {
          stop(paste0("No imputation model possible as a part of the observations that have missing values in block", current_block, " have no other values in [another block] for that other observations exist that have values for the current block and for [another block]."))
        }# else {
        #   # restrict the x values/observations used for the imputation model
        #   # to the blocks for which the current missing values have information
        #   x_values <- x_values[, unlist(blocks[blocks_without_missing])]
        # }
        
        
        # first, take all observations which have values for the maximum number
        # of blocks. If this is below treshold.available.cases, drop the last
        # block and calculate the number of observations again. If the number is
        # still below the threshold, drop the next block etc. until the block
        # with the highest importance is left (smallest number/most left)
        # if this still doesn't satisfy the threshold, drop the block with the
        # highest importance and repeat the steps above
        
        # get all possibilities how the blocks can be combined if there is more
        # than one block
        if (length(blocks_without_missing) > 1) {
          true_false_per_block <- lapply(seq_along(blocks_without_missing),
                                         function(i) {
                                           c(TRUE, FALSE)
                                         })
          block_combinations <- expand.grid(true_false_per_block)
          block_combinations <- block_combinations[, ncol(block_combinations):1]
          # the last combination doesn't contain any block, drop it
          block_combinations <- block_combinations[1:(nrow(block_combinations) - 1), ]
          
          # convert TRUE/FALSE into the actual block numbers
          block_combinations <- apply(block_combinations, 1, function(i) {
            blocks_without_missing[as.logical(i)]
          })
        } else {
          # there is only one possible block
          block_combinations <- blocks_without_missing
        }
        
        # try out all combinations of blocks
        n_obs_information <- data.frame(index = 1:length(block_combinations),
                                        n_obs = rep(0, length(block_combinations)),
                                        greater_threshold = rep(FALSE, length(block_combinations)))
        for (i in 1:length(block_combinations)) {
          possible_cases_index <-
            complete.cases(x_values[, unlist(blocks[block_combinations[[i]]])])
          n_obs_information$n_obs[i] <- sum(possible_cases_index)
          n_obs_information$greater_threshold[i] <- sum(possible_cases_index) >= 
            mcontrol$threshold.available.cases
        }
        
        if (sum(n_obs_information$greater_threshold) == 0) {
          stop(paste0("For block ", current_block, ", no imputation model could be fitted as the number of observations is lower than threshold.available.cases = ", mcontrol$threshold.available.cases))
        }
        
        if (mcontrol$select.available.cases == "maximise.blocks") {
          # only take the combinations which are >= threshold
          temp <- n_obs_information[n_obs_information$greater_threshold, ]
          # the first row is the observation with the best combination of blocks
          index_used_blocks <- temp$index[1]
        }
        if (mcontrol$select.available.cases == "max") {
          # take the combination that leads to the max. number of observations
          index_used_blocks <- which.max(n_obs_information$n_obs)
        }
        # get the corresponding x values
        possible_cases_index <-
          complete.cases(x_values[, unlist(blocks[block_combinations[[index_used_blocks]]])])
        # ncol(x_values) choses the column row_index which was added as the last
        # column further up
        x_values <- x_values[possible_cases_index,
                             c(unlist(blocks[block_combinations[[index_used_blocks]]]),
                               ncol(x_values))]
        blocks_used_for_imputation <- block_combinations[[index_used_blocks]]
        message(paste0("For block ", current_block, " the imputation model contains ", n_obs_information$n_obs[index_used_blocks], " observations."))
        
        
      }
      
      # get the y values (offsets) for the imputation
      # for complete.cases:
      # only take the offsets from the observations that are complete cases
      # over all blocks
      # for available.cases:
      # only take the offsets from the observations that are used for the
      # imputation model (for which the x values are taken)
      y_values_index <- as.vector(calculated_offsets[, 2]) %in%
        as.vector(x_values[, ncol(x_values)])
      y_values_index_num <- calculated_offsets[y_values_index, 2]
      y_values <- calculated_offsets[y_values_index, 1]
      
      # check that the indices of the x and y values are the same
      if (!all.equal(y_values_index_num, x_values[, ncol(x_values)])) {
        stop("Mismatch of covariates and offsets in the imputation model.")
      }
      # delete the index column
      x_values <- x_values[, -ncol(x_values)]
      
      # if the offsets (y_values) are all equal, an imputation
      # model is not possible; instead use the constant offset also for the rest
      # of the imputations
      results <- tryCatch({
        imputation_model <- cv.glmnet(x = x_values,
                                      y = y_values,
                                      nfolds = mcontrol$nfolds.imputation)
        if (!is.null(current_missings)) {
          if (mcontrol$impute.offset.cases == "complete.cases") {
            new_x <- X[current_missings,
                       -blocks[[current_block]]]
          }
          if (mcontrol$impute.offset.cases == "available.cases") {
           new_x <- X[current_missings,
                      unlist(blocks[blocks_used_for_imputation])] 
          }
          
          missing_offsets <- predict(imputation_model,
                                     newx = new_x,
                                     s = mcontrol$lambda.imputation)
          ret <- list(imputation_model = imputation_model,
                      missing_offsets = missing_offsets)
        } else {
          ret <- list(imputation_model = imputation_model)
        }
        ret
        
      }, error = function(e) {
        error_is_constant_issue <- grepl(
          pattern = "y is constant; gaussian glmnet fails at standardization step",
          x = e)
        if (error_is_constant_issue) {
          warning("The offsets calculated for the current block are all equal. An imputation model is not possible, instead the value of the calculated offsets is used as the imputed value.")
          imputation_model <- mean(y_values)
          if (!is.null(current_missings)) {
            missing_offsets <- rep(mean(y_values),
                                   times = length(current_missings))
            ret <- list(imputation_model = imputation_model,
                        missing_offsets = missing_offsets)
          } else {
            ret <- list(imputation_model = imputation_model)
          }
        } else {
          warning(paste0("An error in the imputation model for block ",
                         current_block, " occured."))
          imputation_model <- NULL
          if (!is.null(current_missings)) {
            missing_offsets <- rep(NA, times = length(current_missings))
            ret <- list(imputation_model = imputation_model,
                        missing_offsets = missing_offsets)
          } else {
            ret <- list(imputation_model = imputation_model)
          }
        }
        ret
      })
      
      imputation_model <- results[["imputation_model"]]
      
      # add the observation index
      if (!is.null(missing_offsets)) {
        missing_offsets <- cbind(results[["missing_offsets"]], current_missings)
      }
    }
    if (is.null(current_missings)) {
      new_offsets <- calculated_offsets
    } else {
      new_offsets <- rbind(calculated_offsets, missing_offsets)
    }
    # bring everything into the correct order
    index_sorting <- order(new_offsets[, 2])
    new_offsets <- new_offsets[index_sorting, 1]
  }
  
  # return the complete offsets, the imputation model and the blocks used for
  # the imputation model
  list(new_offsets = new_offsets,
       imputation_model = imputation_model,
       blocks_used_for_imputation = blocks_used_for_imputation)
}
