#' Calculates the offsets for the current block
#'
#' @param current_missings index vector (indices) of current missing observations
#' @param current_observations index vector (indices) of current used observations
#' @param mcontrol control for missing data handling
#' @param current_block index of current block
#' @param pred predictions of current block
#' @param liste list with offsets
#'
#' @return vector of offsets
calculate_offsets <- function(current_missings,
                              current_observations,
                              mcontrol,
                              current_block,
                              pred,
                              liste) {
  
  # store the results for the current block
  # calculate the new offsets
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
      # TODO: do stuff here
    }
    new_offsets <- rbind(calculated_offsets, old_offsets)
    # bring everything into the correct order
    index_sorting <- order(new_offsets[, 2])
    new_offsets <- new_offsets[index_sorting, 1]
  }
  
  # return the vector of the complete offsets
  new_offsets
}