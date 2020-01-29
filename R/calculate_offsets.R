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
#'
#' @return vector of offsets
calculate_offsets <- function(current_missings,
                              current_observations,
                              mcontrol,
                              current_block,
                              pred,
                              liste,
                              X,
                              blocks) {
  
  # store the results for the current block
  # calculate the new offsets
  browser()
  if (is.null(current_missings)) {
    new_offsets <- as.matrix(pred)
  } else {
    # the calculated offsets for this block where there are observations
    calculated_offsets <- as.matrix(pred)
    calculated_offsets <- cbind(calculated_offsets, current_observations)
    # if chosen, ignore the missing block
    if (mcontrol$handle.missingdata == "ignore") {
      # for the missing values, take the offset from the previous block
      # if the missings are in the first block, use 0 as offset
      #JH -> needs further clarifaction
      if (current_block == 1) {
        old_offsets <- cbind(rep(0, length(current_missings)),
                             current_missings)
      } else {
        old_offsets <- cbind(liste[[current_block]][current_missings],
                             current_missings)
      }
    }
    # if chosen impute the missing offsets
    if (mcontrol$handle.missingdata == "impute") {
      if (mcontrol$impute.offset.cases == "complete.cases") {
        # get the x values for the imputation
        x_values <- cbind(X, 1:nrow(X))
        # only take the complete cases over all blocks
        x_values <- x_values[complete.cases(x_values), ]
        # for this imputation model, exclude the information about the current
        # block (because this would completeley determine the offset)
        x_values <- x_values[, -blocks[[current_block]]]
        
        # get the y values (offsets) for the imputation
        # only take the offsets from the observations that are complete cases
        # over all blocks
        y_values_index <- as.vector(calculated_offsets[, 2]) %in%
          as.vector(x_values[, ncol(x_values)])
        y_values <- calculated_offsets[y_values_index, 1]
        
        # check that the indices of the x and y values are the same
        if (!all.equal(y_values_index, x_values[, ncol(x_values)])) {
          stop("Mismatch of covariates and offsets in imputation model for complete cases")
        }
        # delete the index column
        x_values <- x_values[, -ncol(x_values)]
      }
      
      # perform the imputation
      imputation_model <- cv.glmnet(x = x_values,
                                    y = y_values,
                                    nfolds = mcontrol$nfolds.imputation)
    }
    new_offsets <- rbind(calculated_offsets, old_offsets)
    # bring everything into the correct order
    index_sorting <- order(new_offsets[, 2])
    new_offsets <- new_offsets[index_sorting, 1]
  }
  
  # return the vector of the complete offsets
  new_offsets
}
