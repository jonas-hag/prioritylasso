coef.prioritylasso <- function(model) {
  if (!"prioritylasso" %in% class(model)) {
    stop("The model has to be of class prioritylasso.")
  }
  
  # get the coefficients
  if (model$call$family == "cox") {
    coefficients <- model$coefficients
  } else {
    if (is.null(model$block1unpen)) {
      coefficients <- model$coefficients
    } else {
      coefficients <- model$coefficients[-1]
    }
  }
  
  # get the intercepts
  if (model$call$family == "cox") {
    intercepts <- NULL
  } else {
    if (is.null(model$block1unpen)) {
      intercepts <- c()
      for (i in 1:length(model$glmnet.fit)) {
        intercepts[i] <- model$glmnet.fit[[i]]$a0[model$lambda.ind[[i]]]
      }
    } else {
      intercepts <- model$coefficients[1]
      for (i in 2:length(model$glmnet.fit)) {
        intercepts[i] <- model$glmnet.fit[[i]]$a0[model$lambda.ind[[i]]]
      }
    }
  }
  
  list(coefficients = coefficients,
       intercepts = intercepts)
}