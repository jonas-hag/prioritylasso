make_imputation_model <- function(x_values,
                                  y_values,
                                  mcontrol,
                                  X,
                                  current_missings,
                                  blocks,
                                  current_block,
                                  blocks_used_for_imputation) {
  results <- tryCatch({
    imputation_model <- cv.glmnet(x = x_values,
                                  y = y_values,
                                  nfolds = mcontrol$nfolds.imputation)
    if (!is.null(current_missings)) {
      if (mcontrol$impute.offset.cases == "complete.cases") {
        new_x <- X[current_missings,
                   -blocks[[current_block]], drop = FALSE]
      }
      if (mcontrol$impute.offset.cases == "available.cases") {
        new_x <- X[current_missings,
                   unlist(blocks[blocks_used_for_imputation]), drop = FALSE] 
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
    # if the offsets (y_values) are all equal, an imputation
    # model is not possible; instead use the constant offset also for the rest
    # of the imputations
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
                     current_block, " occured:"))
      print(e)
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
  
  results
}