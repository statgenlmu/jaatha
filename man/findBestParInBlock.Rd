\name{findBestParInBlock}
\alias{findBestParInBlock}
\title{Functions that uses the fitted GLMs to estimate the parameters that have the
highest likelihood.}
\usage{
findBestParInBlock(block, glm.fitted, sum.stats)
}
\arguments{
  \item{block}{The block we are in.}

  \item{glm.fitted}{The fitted GlMs.}

  \item{sum.stats}{The observed summary statistics}
}
\value{
A list with maximum likelihood parameter (est) and
log-likelihood (score)
}
\description{
Functions that uses the fitted GLMs to estimate the
parameters that have the highest likelihood.
}

