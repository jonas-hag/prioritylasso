#' Compare the rows of a matrix with a pattern
#'
#' @param object matrix
#' @param pattern pattern which is compared against the rows of the matrix
#'
#' @return logical vector if the pattern matches the rows
#' @importFrom checkmate assert_matrix
compare_boolean <- function(object,
                            pattern) {
  pattern <- as.vector(pattern)
  assert_matrix(object, ncols = length(pattern))
  apply(object, 1, function(x) {
    x <- as.vector(x)
    res <- all.equal(x, pattern)
    if (!is.logical(res)) {
      res <- FALSE
    }
    res
  })
}
