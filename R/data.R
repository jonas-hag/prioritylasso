#' Simulated AML data with binary outcome
#'
#' A data set containing the binary outcome and 1028 predictor variables of 400 artificial AML patients.
#'
#' We generated the data in the following way: We took the empirical correlation of 1028 variables related to
#' 315 AML patients. This correlation served as a correlation matrix when generating 1028 multivariate
#' normally distributed variables with the R function \code{\link[mvtnorm]{rmvnorm}}. Because we didn't have a positive
#' definite matrix, we took the nearest positive definite matrix according to the function \code{\link[Matrix]{nearPD}}.
#' The variables that should be binary were dichotomized, so that their marginal probabilities corresponded to
#' the marginal probabilities they were based on.
#' The coefficients were defined by
#' \itemize{
#' \item{\code{beta_b1 <- c(0.8, 0.8, 0.6, 0.6)}}
#' \item{\code{beta_b2 <- c(rep(0.5,3), rep(0,2))}}
#' \item{\code{beta_b3 <- c(rep(0.4, 4), rep(0,15))}}
#' \item{\code{beta_b4 <- c(rep(0.5, 5), rep(0.3, 5), rep(0,990))}}.
#' }
#' We included them in the vector \code{beta <- c(beta_b1, beta_b2, beta_b3, beta_b4)} and calculated
#' the probability through \deqn{pi = exp(\beta*x) / (1 + exp(\beta*x))} where x denotes our data matrix
#' with 1028 predictor variables. Finally we got the outcome through
#' \code{pl_out <- rbinom(400, size = 1, p = pi)}.
#'
#'
#'
#' @format A data frame with 400 rows and 1029 variables:
#' \describe{
#'   \item{pl_out: (\code{pl_data[,1029]})}{binary outcome representing refractory status.}
#'   \item{b1: (\code{pl_data[,1:4]})}{4 binary variables representing variables with a known influence on the outcome.}
#'   \item{b2: (\code{pl_data[,5:9]})}{5 continuous variables representing clinical variables.}
#'   \item{b3: (\code{pl_data[,10:28]})}{19 binary variables representing mutations.}
#'   \item{b4: (\code{pl_data[,29:1028]})}{1000 continuous variables representing gene expression data.}
#' }
#'
#'
"pl_data"
