\name{Jaatha.confidenceIntervals}
\alias{Jaatha.confidenceIntervals}
\title{Function for calculating bootstrap confidence intervals for
jaatha estimates.}
\usage{
  Jaatha.confidenceIntervals(jaatha, conf.level = 0.95,
    replicas = 100, cores = 1,
    log.folder = paste(tempdir(), "jaatha-logs", sep = "/"))
}
\arguments{
  \item{jaatha}{A jaatha object that was returned by
  Jaatha.refinedSearch.}

  \item{conf.level}{The intended confidence level of the
  interval. The actual level can vary slightly.}

  \item{replicas}{The number of Jaatha runs that we perform
  for calculating the confidence interval. Should be
  reasonable large.}

  \item{cores}{The number of CPU cores that will be used.}

  \item{log.folder}{A folder were log-files for the
  indiviual runs are placed.}
}
\value{
  The Jaatha Object with confidence intervals included.
}
\description{
  This functions calculates bias-corrected and accelerated
  (BCa) bootstrap confidence intervals as described in
  Section 14.3 of "An introduction to the bootstrap" by
  Efron and Tibshirani. Basically, we simulate many
  datasets under the model using the estimated parameters
  and do a complete Jaatha estimation procedure on each
  dataset. We can then use the thereby estimated values to
  calculate approximate confidence intervals.
}
\details{
  Warning: This requires a large number of Jaatha runs and
  should be best executed on a small cluster rather than on
  a desktop computer.
}

