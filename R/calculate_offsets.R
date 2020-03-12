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
  missingness_pattern <- NULL
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
        # store the missingness pattern, TRUE means missing
        missingness_pattern <- seq_along(blocks) %in% current_block
        
        # get the y values (offsets) for the imputation
        # for complete.cases:
        # only take the offsets from the observations that are complete cases
        # over all blocks
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
        
        # generate the imputation model (if the y_values are constant, the mean
        # of these values is used)
        results <- make_imputation_model(x_values = x_values,
                                         y_values = y_values,
                                         mcontrol = mcontrol,
                                         X = X,
                                         current_missings = current_missings,
                                         blocks = blocks,
                                         current_block = current_block,
                                         blocks_used_for_imputation = 
                                           blocks_used_for_imputation)
        
        imputation_model <- results[["imputation_model"]]
        if (!is.null(current_missings)) {
          # add the observation index
          missing_offsets <- cbind(results[["missing_offsets"]], current_missings)
        }
      }
      
      if (mcontrol$impute.offset.cases == "available.cases" &&
          !is.null(current_missings)) {
        # get the x values for the imputation
        x_values <- cbind(X, row_index = 1:nrow(X))
        # only take the observations which have values for the current block
        # and other blocks
        
        # get all the x values that have information for the current block
        row_index <- complete.cases(x_values[, blocks[[current_block]]])
        x_values <- x_values[row_index, ]
        
        # get the information which other blocks of the observations which are
        # missing the current block are missing
        missing_index_overview <- matrix(FALSE, nrow = length(current_missings),
                                         ncol = length(blocks))
        for (i in seq_along(blocks)) {
          missing_index_overview[, i] <- !complete.cases(X[current_missings, blocks[[i]]])
        }
        
        # determine the unique pattern of missingness
        unique_missingness_pattern <- unique(missing_index_overview)
        # determine an imputation model for every unique missingness pattern
        results_available_cases <-
          apply(unique_missingness_pattern, 1, function(pattern) {
            # determine all observerations (in which the current block is
            # missing) that have the same missingness pattern
            index_missingness_pattern <- compare_boolean(missing_index_overview,
                                                         pattern)
            missings_current_pattern <- current_missings[index_missingness_pattern]
            
            # calculate an imputation model for the observations with this
            # missingness pattern
            results <- make_imputation_available_cases(X = X,
                                                       current_missings = missings_current_pattern,
                                                       blocks = blocks,
                                                       current_block = current_block,
                                                       mcontrol = mcontrol,
                                                       x_values = x_values,
                                                       calculated_offsets = calculated_offsets,
                                                       pattern_miss = pattern)
            list(results = results,
                 used_missing_observations = missings_current_pattern,
                 missingness_pattern = pattern)
            
          })
        
      imputation_model <- lapply(results_available_cases, function(x) {
        x$results[["imputation_model"]]
      })
      missingness_pattern <- lapply(results_available_cases, function(x) {
        x$missingness_pattern
      })
      blocks_used_for_imputation <- lapply(results_available_cases, function(x) {
        x$results[["blocks_used_for_imputation"]]
      })
      if (!is.null(current_missings)) {
        # add the observation index
        missing_offsets <- lapply(results_available_cases, function(x) {
          cbind(x$results[["missing_offsets"]], x$used_missing_observations)
        })
        missing_offsets <- do.call("rbind", missing_offsets)
      }
        
        
      }
      
    }
    if (is.null(current_missings)) {
      new_offsets <- calculated_offsets
    } else {
      # bring the already calculated offsets and the imputed offsets together
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
       blocks_used_for_imputation = blocks_used_for_imputation,
       missingness_pattern = missingness_pattern)
}
