make_imputation_available_cases <- function(X,
                                            current_missings,
                                            blocks,
                                            current_block,
                                            mcontrol,
                                            x_values,
                                            calculated_offsets,
                                            pattern_miss) {
  # get further information which blocks can be used for the imputation
  # model
  # -> get all the blocks for which the current missing values have
  # information
  na_matrix_current_missings <- is.na(X[current_missings, , drop = FALSE])
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
  
  # this safety condition shouldn't be needed anymore
  if (length(blocks_without_missing) == 0) {
    stop(paste0("For block ", current_block, " and missingness pattern ",
                paste0(pattern_miss, collapse = " "), "no imputation model is possible, as there is no other block where data is not missing."))
  }
  
  
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
  message(paste0("For block ", current_block, " and the missingness pattern ",
                 paste0(pattern_miss, collapse = " ")," the imputation model contains ",
                 n_obs_information$n_obs[index_used_blocks], " observations."))
  
  # get the y values (offsets) for the imputation
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
    missing_offsets <- results[["missing_offsets"]]
  } else {
    missing_offsets <- NULL
  }
  
  # return the results
  list(imputation_model = imputation_model,
       missing_offsets = missing_offsets,
       blocks_used_for_imputation = blocks_used_for_imputation)
}