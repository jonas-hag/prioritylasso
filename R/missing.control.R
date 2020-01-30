#' Construct control structures for handling of missing data for \code{prioritylasso}
#'
#' @param handle.missingdata how blockwise missing data should be treated. Default is \code{none} which does nothing, \code{ignore} ignores the observations with missing data for the current block, \code{impute.offset} imputes the offset for the missing values.
#' @param impute.offset.cases which cases should be used for the imputation model to impute missing offsets. So far, only complete cases are supported
#' @param nfolds.imputation nfolds for the glmnet of the imputation model
#' @param lambda.imputation which lambda-value should be used for predicting the imputed offsets in cv.glmnet
#' @param perc.comp.cases.warning percentage of complete cases when a warning is issued of too few cases for the imputation model
#'
#' @return list with control parameters
#' @export
#' @importFrom checkmate assert_numeric assert_number
missing.control <- function(handle.missingdata = c("none", "ignore",
                                                   "impute.offset"),
                            impute.offset.cases = "complete.cases",
                            nfolds.imputation = 10,
                            lambda.imputation = c("lambda.min", "lambda.1se"),
                            perc.comp.cases.warning = 0.3) {
  
  handle.missingdata <- match.arg(handle.missingdata)
  impute.offset.cases <- match.arg(impute.offset.cases)
  assert_number(nfolds.imputation, lower = 3)
  lambda.imputation <- match.arg(lambda.imputation)
  assert_numeric(perc.comp.cases.warning, lower = 0, upper = 1, max.len = 1)
  
  result <- list(handle.missingdata = handle.missingdata,
                 impute.offset.cases = impute.offset.cases,
                 nfolds.imputation = nfolds.imputation,
                 lambda.imputation = lambda.imputation,
                 perc.comp.cases.warning = perc.comp.cases.warning)
  
  class(result) <- c("pl.missing.control", class(result))
  result
}
