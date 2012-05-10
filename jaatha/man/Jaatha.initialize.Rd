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

  \item{resultFile}{ABC}

  \item{debugMode}{def}

  \item{logFile}{geh}
}
\value{
  A S4-Object of type jaatha containing the settings
}
\description{
  This function sets the basic parameters for an analysis
  with Jaatha and is the first step for each application of
  it.
}

