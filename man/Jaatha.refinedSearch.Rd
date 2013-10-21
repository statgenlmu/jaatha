\name{Jaatha.refinedSearch}
\alias{Jaatha.refinedSearch}
\title{Iterative search for the maximum composite likelihood parameters}
\usage{
  Jaatha.refinedSearch(jaatha, best.start.pos, sim,
    sim.final, epsilon = 0.2, half.block.size = 0.05,
    weight = 0.9, max.steps = 200)
}
\arguments{
  \item{jaatha}{The Jaatha settings (create with
  \code{\link{Jaatha.initialize}})}

  \item{best.start.pos}{This is the number of best starting
  positions found in the initial search that we will use.
  Jaatha runs a separate search starting from each of this
  points.}

  \item{sim}{The number of simulations that are performed
  in each step}

  \item{sim.final}{The number of simulations that are
  performed after the search to estimate the composite log
  likelihood. If not specified, the value of \code{sim}
  will be used}

  \item{epsilon}{The search stops if the improvement of the
  score is less than this for 5 times in a row.}

  \item{half.block.size}{The size of the new block that is
  created around a new maximum.}

  \item{weight}{The weighting factor that will reduce the
  influence of old block in the estimation procedure}

  \item{max.steps}{The search will stop at this number of
  steps if not stopped before (see \code{epsilon})}
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

