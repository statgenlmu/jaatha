\name{Jaatha.initialize}
\alias{Jaatha.initialize}
\title{Basic initialization of a Jaatha estimation}
\usage{
  Jaatha.initialize(demographicModel, sumStats,
    seed = numeric(), resultFile = "", debugMode = F,
    logFile = "")
}
\arguments{
  \item{demographicModel}{The demographic model to use}

  \item{sumStats}{The summary statistics calculated from
  the real data}

  \item{seed}{An integer used as seed for both Jaatha and
  the simulation software}

  \item{resultFile}{A File in which the results of the
  search will be written.}

  \item{debugMode}{If TRUE, a debug output will be
  produced}

  \item{logFile}{If specified, the debug output will be
  redirected to this file}
}
\value{
  A S4-Object of type jaatha containing the settings
}
\description{
  This function sets the basic parameters for an analysis
  with Jaatha and is the first step for each application of
  it.
}

