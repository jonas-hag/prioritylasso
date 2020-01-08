#' Predictions from prioritylasso
#'
#' Makes predictions for a \code{prioritylasso} object. It can be chosen between linear predictors or fitted values.
#'
#' @param object An object of class \code{prioritylasso}.
#' @param newdata (nnew \code{x} p) matrix or data frame with new values.
#' @param type Specifies the type of predictions. \code{link} gives the linear predictors for all types of response and \code{response} gives the fitted values.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Predictions that depend on \code{type}.
#'
#' @author Simon Klau
#' @export
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
                                  ...){

  if (is.null(newdata)) {
    newdata <- as.matrix(object$X)
  } else {
    newdata <- data.matrix(newdata)
  }

  type <- match.arg(type)

  if(is.null(object$block1unpen)){
    if(object$call$family == "cox"){
      pred <- newdata %*% object$coefficients
    } else {
      pred <- newdata %*% object$coefficients + object$glmnet.fit[[1]]$a0[object$lambda.ind[[1]]]
    }
  } else {
    if(object$call$family == "cox"){
      pred <- newdata %*% object$coefficients
    } else {
      coeff <- object$coefficients[-1]
      intercept <- object$coefficients[1]

      pred <- newdata %*% coeff + intercept
    }
  }


  if(type == "response"){
    if(object$call$family == "binomial"){
      pred <- exp(pred)/(1+exp(pred)) # fitted probabilities
    }
    if(object$call$family == "cox"){
      pred <- exp(pred) # fitted relative risk (risk score (exp(lp)))
    }
  }
  return(pred)
}
