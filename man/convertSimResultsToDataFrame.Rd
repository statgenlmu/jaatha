\name{convertSimResultsToDataFrame}
\alias{convertSimResultsToDataFrame}
\title{Converts simulation results into a data frame that is usable for fitting a
glm.}
\usage{
convertSimResultsToDataFrame(sim.data, sum.stat)
}
\arguments{
  \item{sim.data}{Results from simulations}

  \item{sum.stat}{Name of the summary statistics which
  should get converted}
}
\value{
The summary statistics as data.frame
}
\description{
Currently only works with nx2 matix summary statistics and
vectorizes thoose.
}

