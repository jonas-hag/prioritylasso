# prioritylasso

<!-- badges: start -->
[![R-CMD-check](https://github.com/jonas-hag/prioritylasso/workflows/R-CMD-check/badge.svg)](https://github.com/jonas-hag/prioritylasso/actions)
[![R-CMD-check](https://github.com/jonas-hag/prioritylasso/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jonas-hag/prioritylasso/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

prioritylasso is an `R` package that fits successive Lasso models for several
blocks of (omics) data with different priorities and takes the predicted values
as an offset for the next block.

## Installation
For the latest stable release from CRAN use:

```r
install.packages("prioritylasso")
```

To get the latest version from github, use:

```r
remotes::install_github("jonas-hag/prioritylasso")
```

## Usage
The basic functionality is provided by the `prioritylasso` function. You can run
a simple model with a gaussian dependent variable:

```r
results <- prioritylasso(
  X = matrix(rnorm(50 * 500), 50, 500),
  Y = rnorm(50),
  family = "gaussian",
  type.measure = "mse",
  blocks = list(bp1 = 1:75, bp2 = 76:200, bp3 = 201:500),
  max.coef = c(Inf, 8, 5),
  block1.penalization = TRUE,
  lambda.type = "lambda.min",
  standardize = TRUE,
  nfolds = 5,
  cvoffset = FALSE
)
```

Binary outcome data and Cox models are also possible. For a better overview,
have a look at the [introductory vignette](vignettes/prioritylasso_vignette.Rmd).

### Block-wise missing data
A special type of missing data is block-wise missing data and occurs when the
data contains "blocks", e.g. several variables that belong together like
clinical measurements, mRNA sequencing data, SNP data etc. This means that for
some observations not all blocks are observed. To deal with this type of
missingness, prioritylasso provides the following options to fit a model to a
data set:

- `ignore`: the Lasso model for every block is only fitted
with the observations that have no missing values for this block. For
observations with the current block missing, the offset from the previous
block is carried forward
- `impute`: the Lasso model for every block is only fitted
with the observations that have no missing values for this block. For
observations with the current block missing, the offset from the previous
block is imputed. The imputation model is either based on all other blocks or
it is tried to use as much information as possible for more complex missingness
patterns.

These options can be set in the function `missing.control`.

If a prioritylasso model should be used to predict on data with block-wise missing
data, the following options are available:

- `set.zero`: ignores the missing data for the calculation of the prediction
(the missing value is set to zero)
- `impute.block`: use an imputation model to impute the offset of a missing
block. In order to work, the prioritylasso model must be trained with the option
`impute` and the missingness patterns in the test data have to be the same as in
the train data

These options can be set in `handle.missingtestdata` of the `predict` function.

## Paper
For more information about the method, see the following paper:

>[Klau, Simon, Vindi Jurinovic, Roman Hornung, Tobias Herold, and Anne-Laure Boulesteix. "Priority-Lasso: a simple hierarchical approach to the prediction of clinical outcome using multi-omics data." BMC bioinformatics 19, no. 1 (2018): 1-14.](https://doi.org/10.1186/s12859-018-2344-6)
