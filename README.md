# prioritylasso
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

## Paper
For more information about the method, see the following paper:

>[Klau, Simon, Vindi Jurinovic, Roman Hornung, Tobias Herold, and Anne-Laure Boulesteix. "Priority-Lasso: a simple hierarchical approach to the prediction of clinical outcome using multi-omics data." BMC bioinformatics 19, no. 1 (2018): 1-14.](https://doi.org/10.1186/s12859-018-2344-6)
