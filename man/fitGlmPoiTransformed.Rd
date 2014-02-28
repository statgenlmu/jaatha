\name{fitGlmPoiTransformed}
\alias{fitGlmPoiTransformed}
\title{Fits a GLM for a summary statistics of type "poisson.transformed" and
"poisson.independent"}
\usage{
fitGlmPoiTransformed(sim.data, sum.stat, transformation, weighting, jaatha)
}
\arguments{
  \item{sim.data}{Results from simulations}

  \item{sum.stat}{Name of the summary statistics}

  \item{transformation}{Transformation function to apply to
  the data}

  \item{weighting}{Potentially weights for the
  simulations.}

  \item{jaatha}{A Jaatha Object.}
}
\value{
A list of fitted GLMs, one for each function
}
\description{
Fits a GLM for a summary statistics of type
"poisson.transformed" and "poisson.independent"
}

