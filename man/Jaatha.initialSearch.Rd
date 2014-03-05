\name{Jaatha.initialSearch}
\alias{Jaatha.initialSearch}
\title{Search the parameter space for good starting positions}
\usage{
Jaatha.initialSearch(jaatha, sim = 200, blocks.per.par = 2, rerun = FALSE)
}
\arguments{
  \item{jaatha}{The Jaatha settings (create with
  \code{\link{Jaatha.initialize}})}

  \item{sim}{Numeric. The number of simulations that are
  performed in each bin}

  \item{blocks.per.par}{Numeric. The number of block per
  parameter.  Will result in Par^blocks.per.par blocks}

  \item{rerun}{You can repeat a previously done initial
  search in Jaatha.  Do do so, just call the initial search
  function with the jaatha object result of the first
  initial search and set rerun to 'TRUE'.}
}
\value{
The jaatha object with starting positions
}
\description{
This functions devides the parameter space in different
parts (blocks). In each block, simulations for different
parameter combinations are run to roughly predict the
combination with the highest score (which is equivalent to
the highest composite log likelihood). This points can
later be used as starting positions for the secound
estimation phase of Jaatha
(\code{\link{Jaatha.refinedSearch}})
}

