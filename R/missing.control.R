#' Construct control structures for handling of missing data for \code{prioritylasso}
#'
#' @param handle.missingdata how blockwise missing data should be treated. Default is \code{none} which does nothing, \code{ignore} ignores the observations with missing data for the current block, \code{impute.offset} imputes the offset for the missing values.
#' @param impute.offset.cases which cases should be used for the imputation model to impute missing offsets. So far, only complete cases are supported
#' @param perc.comp.cases.warning percentage of complete cases when a warning is issued of too few cases for the imputation model
#'
#' @return list with control parameters
#' @export
#' @importFrom checkmate assert_numeric
missing.control <- function(handle.missingdata = c("none", "ignore",
                                                   "impute.offset"),
                            impute.offset.cases = "complete.cases",
                            perc.comp.cases.warning = 0.1) {
  
  handle.missingdata <- match.arg(handle.missingdata)
  impute.offset.cases <- match.arg(impute.offset.cases)
  assert_numeric(perc.comp.cases.warning, lower = 0, upper = 1, max.len = 1)
  
  result <- list(handle.missingdata = handle.missingdata,
                 impute.offset.cases = impute.offset.cases,
                 perc.comp.cases.warning = perc.comp.cases.warning)
  
  class(result) <- c("pl.missing.control", class(result))
  result
}