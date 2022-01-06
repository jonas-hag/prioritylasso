#' Extract coefficients from a prioritylasso object
#'
#' @param object model of type prioritylasso
#' @param ... additional arguments, currently not used
#'
#' @return List with the coefficients and the intercepts
#' @export
coef.prioritylasso <- function(object, ...) {
  if (!"prioritylasso" %in% class(object)) {
    stop("The model has to be of class prioritylasso.")
  }
  
  # get the coefficients
  if (object$call$family == "cox") {
    coefficients <- object$coefficients
  } else {
    if (is.null(object$block1unpen)) {
      coefficients <- object$coefficients
    } else {
      coefficients <- object$coefficients[-1]
    }
  }
  
  # get the intercepts
  if (object$call$family == "cox") {
    intercepts <- NULL
  } else {
    if (is.null(object$block1unpen)) {
      intercepts <- c()
      for (i in 1:length(object$glmnet.fit)) {
        intercepts[i] <- object$glmnet.fit[[i]]$a0[object$lambda.ind[[i]]]
      }
    } else {
      intercepts <- object$coefficients[1]
      for (i in 2:length(object$glmnet.fit)) {
        intercepts[i] <- object$glmnet.fit[[i]]$a0[object$lambda.ind[[i]]]
      }
    }
  }
  
  list(coefficients = coefficients,
       intercepts = intercepts)
}