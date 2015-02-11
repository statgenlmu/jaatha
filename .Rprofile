.First <- function() {
  .First.sys()                            # Load the default packages
  devtools::load_all(".", quiet = TRUE)   # Source the jaatha package
  if(!require('testthat', quietly = TRUE)) warning('testthat not installed')
  source('tests/testthat/setup-custom-simulation-interface.R')
}

