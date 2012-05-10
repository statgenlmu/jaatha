\name{Jaatha.refineSearch}
\alias{Jaatha.refineSearch}
\title{Iterative search for the maximum composite likelihood parameters}
\usage{
  Jaatha.refineSearch(jObject, startPoints, nSim,
    nFinalSim, epsilon = 0.2, halfBlockSize = 0.05,
    weight = 0.9, nMaxStep = 200)
}
\arguments{
  \item{jObject}{The Jaatha settings (create with
  \code{\link{Jaatha.initialize}})}

  \item{startPoints}{From each of starting positions in
  this list, a refined search will be performed starting in
  this position.  (create this list with
  \code{\link{Jaatha.initialSearch}})}

  \item{nSim}{The number of simulations that are performed
  in each step}

  \item{nFinalSim}{The number of simulations that are
  performed after the search to estimate the composite log
  likelihood. If not specified, the value of \code{nSim}
  will be used}

  \item{epsilon}{The search stops if the improvement of the
  score is less than this for 5 times in a row.}

  \item{halfBlockSize}{The size of the new block that is
  created around a new maximum.}

  \item{weight}{The weighting factor that will reduce the
  influence of old block in the estimation procedure}

  \item{nMaxStep}{The search will stop at this number of
  steps if not stopped before (see epsilon)}
}
\value{
  An Jaatha object. The found values are written to the
  slot likelihood.table.
}
\description{
  This function searches for the parameter combination with
  the highest composite likelihood. Therefore it
  iteratively searches the area (=block) around the last
  found value for a new maximum using simulations and a
  generalized linear model.
}

