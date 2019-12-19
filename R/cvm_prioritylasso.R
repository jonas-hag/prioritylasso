#' prioritylasso with several block specifications
#'
#' Runs prioritylasso for a list of block specifications and gives the best results
#' in terms of cv error.
#'
#' @param X a (nxp) matrix or data frame of predictors with observations in rows and predictors in columns.
#' @param Y n-vector giving the value of the response (either continuous, numeric-binary 0/1, or \code{Surv} object).
#' @param weights observation weights. Default is 1 for each observation.
#' @param family should be "gaussian" for continuous \code{Y}, "binomial" for binary \code{Y}, "cox" for \code{Y} of type \code{Surv}.
#' @param type.measure The accuracy/error measure computed in cross-validation. It should be "class" (classification error) or "auc" (area under the ROC curve) if \code{family="binomial"}, "mse" (mean squared error) if \code{family="gaussian"} and "deviance" if \code{family="cox"} which uses the partial-likelihood.
#' @param blocks.list list of the format \code{list(list(bp1=...,bp2=...,), list(bp1=,...,bp2=...,), ...)}. For the specification of the entries, see \code{\link[prioritylasso]{prioritylasso}}.
#' @param max.coef.list list of \code{max.coef} vectors. The first entries are omitted if \code{block1.penalization = FALSE}. Default is \code{NULL}.
#' @param block1.penalization whether the first block should be penalized. Default is TRUE.
#' @param lambda.type specifies the value of lambda used for the predictions. \code{lambda.min} gives lambda with minimum cross-validated errors. \code{lambda.1se} gives the largest value of lambda such that error is within 1 standard error of the minimum. Note that \code{lambda.1se} can only be chosen without restrictions of \code{max.coef}.
#' @param standardize logical, whether the predictors should be standardized or not. Default is TRUE.
#' @param nfolds the number of CV procedure folds.
#' @param foldid an optional vector of values between 1 and nfold identifying what fold each observation is in.
#' @param cvoffset logical, whether CV should be used to estimate the offsets. Default is FALSE.
#' @param cvoffsetnfolds the number of folds in the CV procedure that is performed to estimate the offsets. Default is 10. Only relevant if \code{cvoffset=TRUE}.

#' @param ... Other arguments that can be passed to the function \code{cv.glmnet}.
#'
#' @return object of class \code{prioritylasso} with the following elements. If these elements are lists, they contain the results for each penalized block of the best result.
#' \describe{
#' \item{\code{lambda.ind}}{list with indices of lambda for \code{lambda.type}.}
#' \item{\code{lambda.type}}{type of lambda which is used for the predictions.}
#' \item{\code{lambda.min}}{list with values of lambda for \code{lambda.type}.}
#' \item{\code{min.cvm}}{list with the mean cross-validated errors for \code{lambda.type}.}
#' \item{\code{nzero}}{list with numbers of non-zero coefficients for \code{lambda.type}.}
#' \item{\code{glmnet.fit}}{list of fitted \code{glmnet} objects.}
#' \item{\code{name}}{a text string indicating type of measure.}
#' \item{\code{block1unpen}}{if \code{block1.penalization = FALSE}, the results of either the fitted \code{glm} or \code{coxph} object.}
#' \item{\code{best.blocks}}{character vector with the indices of the best block specification.}
#' \item{\code{best.max.coef}}{vector with the number of maximal coefficients corresponding to \code{best.blocks}.}
#' \item{\code{coefficients}}{coefficients according to the results obtained with \code{best.blocks}.}
#' \item{\code{call}}{the function call.}
#' }
#'
#' @note The function description and the first example are based on the R package \code{ipflasso}.
#' @author Simon Klau \cr
#'         Maintainer: Simon Klau (\email{simonklau@ibe.med.uni-muenchen.de})
#' @references Klau, S., Jurinovic, V., Hornung, R., Herold, T., Boulesteix, A.-L. (2018). Priority-Lasso: a simple hierarchical approach to the prediction of clinical outcome using multi-omics data. BMC Bioinformatics 19, 322
#' @export
#' @seealso \code{\link[prioritylasso]{pl_data}}, \code{\link[prioritylasso]{prioritylasso}}, \code{\link[ipflasso]{cvr2.ipflasso}}
#' @examples
#' cvm_prioritylasso(X = matrix(rnorm(50*500),50,500), Y = rnorm(50), family = "gaussian",
#'                   type.measure = "mse", lambda.type = "lambda.min", nfolds = 5,
#'                   blocks.list = list(list(bp1=1:75, bp2=76:200, bp3=201:500),
#'                                      list(bp1=1:75, bp2=201:500, bp3=76:200)))
#'\dontrun{
#' cvm_prioritylasso(X = pl_data[,1:1028], Y = pl_data[,1029], family = "binomial",
#'                   type.measure = "auc", standardize = FALSE, block1.penalization = FALSE,
#'                   blocks.list = list(list(1:4, 5:9, 10:28, 29:1028),
#'                                      list(1:4, 5:9, 29:1028, 10:28)),
#'                   max.coef.list = list(c(Inf, Inf, Inf, 10), c(Inf, Inf, 10, Inf)))}
#'


cvm_prioritylasso <- function(X, Y, weights, family, type.measure, blocks.list, max.coef.list = NULL,
                              block1.penalization = TRUE, lambda.type = "lambda.min",
                              standardize = TRUE, nfolds = 10, foldid, cvoffset = FALSE, cvoffsetnfolds = 10, ...){

  if(!is.null(max.coef.list)){
    if(length(blocks.list) != length(max.coef.list)){stop("blocks.list and max.coef.list must have the same length.")}
    if(any(unlist(lapply(blocks.list, length)) != unlist(lapply(max.coef.list, length)))){stop("blocks.list and the entries of max.coef.list must have the same length.")}
  }

  nlist <- length(blocks.list)
  all_res <- list()
  cvmin = +Inf

  for(j in 1:nlist){
    all_res[[j]] <- prioritylasso(X, Y, weights, family = family, type.measure = type.measure, blocks = blocks.list[[j]],
                                  max.coef = max.coef.list[[j]], block1.penalization = block1.penalization,
                                  lambda.type = lambda.type, standardize = standardize, nfolds = nfolds, foldid,
                                  cvoffset = cvoffset, cvoffsetnfolds = cvoffsetnfolds, ...)


    if (type.measure != "auc"){
      mincvmj <-  all_res[[j]]$min.cvm[[length(blocks.list[[j]])]]
    }
    if(type.measure == "auc"){
      mincvmj <-  -all_res[[j]]$min.cvm[[length(blocks.list[[j]])]]
    }


    if(mincvmj < cvmin){
      ind.best.pr <- j
      cvmin <- mincvmj
    }
  }


  best.var <- blocks.list[[ind.best.pr]]
  bfv <- lapply(best.var, '[[', 1)
  blv <- lapply(best.var, length)
  best.blocks <- vector()
  for(k in 1:length(bfv)){
    best.blocks <- c(best.blocks, paste("bp",k," = ", bfv[[k]],":",blv[[k]] + bfv[[k]] - 1, sep=""))
  }

  finallist <- list(lambda.ind = all_res[[ind.best.pr]]$lambda.ind, lambda.type = lambda.type,
                    lambda.min = all_res[[ind.best.pr]]$lambda.min,
                    min.cvm = all_res[[ind.best.pr]]$min.cvm, nzero = all_res[[ind.best.pr]]$nzero,
                    glmnet.fit = all_res[[ind.best.pr]]$glmnet.fit, name = all_res[[ind.best.pr]]$name,
                    block1unpen = all_res[[ind.best.pr]]$block1unpen, best.blocks = best.blocks,
                    best.max.coef = max.coef.list[[ind.best.pr]], coefficients = all_res[[ind.best.pr]]$coefficients,
                    call = match.call())

  class(finallist) <- "prioritylasso"

  return(finallist)

}
