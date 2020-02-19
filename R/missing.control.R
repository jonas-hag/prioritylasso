#' Construct control structures for handling of missing data for \code{prioritylasso}
#'
#' @param handle.missingdata how blockwise missing data should be treated. Default is \code{none} which does nothing, \code{ignore} ignores the observations with missing data for the current block, \code{impute.offset} imputes the offset for the missing values.
#' @param offset.firstblock determines if the offset of the first block for missing observations is zero or the intercept of the observed values for \code{handle.missingdata = ignore}
#' @param impute.offset.cases which cases/observations should be used for the imputation model to impute missing offsets. Supported are complete cases (additional constraint is that every observation can only contain one missing block) and all available observations which have an overlap with the current block.
#' @param nfolds.imputation nfolds for the glmnet of the imputation model
#' @param lambda.imputation which lambda-value should be used for predicting the imputed offsets in cv.glmnet
#' @param perc.comp.cases.warning percentage of complete cases when a warning is issued of too few cases for the imputation model
#' @param threshold.available.cases if the number of available cases for \code{impute.offset.cases = available.cases} is below this threshold, \code{prioritylasso} tries to reduce the number of blocks taken into account for the imputation model to increase the number of observations used for the imputation model.
#' @param select.available.cases determines how the blocks which are used for the imputation model are selected when \code{impute.offset.cases = available.cases}. \code{max} selects the blocks that maximise the number of observations, \code{maximise.blocks} tries to include as many blocks as possible, starting with the blocks with the hightes priority 
#'
#' @return list with control parameters
#' @export
#' @importFrom checkmate assert_numeric assert_number
missing.control <- function(handle.missingdata = c("none", "ignore",
                                                   "impute.offset"),
                            offset.firstblock = c("zero", "intercept"),
                            impute.offset.cases = c("complete.cases", "available.cases"),
                            nfolds.imputation = 10,
                            lambda.imputation = c("lambda.min", "lambda.1se"),
                            perc.comp.cases.warning = 0.3,
                            threshold.available.cases = 30,
                            select.available.cases = c("maximise.blocks", "max")) {
  
  handle.missingdata <- match.arg(handle.missingdata)
  offset.firstblock <- match.arg(offset.firstblock)
  impute.offset.cases <- match.arg(impute.offset.cases)
  assert_number(nfolds.imputation, lower = 3)
  lambda.imputation <- match.arg(lambda.imputation)
  assert_number(perc.comp.cases.warning, lower = 0, upper = 1)
  assert_number(threshold.available.cases, lower = 1)
  select.available.cases <- match.arg(select.available.cases)
  
  result <- list(handle.missingdata = handle.missingdata,
                 offset.firstblock = offset.firstblock,
                 impute.offset.cases = impute.offset.cases,
                 nfolds.imputation = nfolds.imputation,
                 lambda.imputation = lambda.imputation,
                 perc.comp.cases.warning = perc.comp.cases.warning,
                 threshold.available.cases = threshold.available.cases,
                 select.available.cases = select.available.cases)
  
  class(result) <- c("pl.missing.control", class(result))
  result
}
