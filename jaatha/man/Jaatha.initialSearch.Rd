\name{Jaatha.initialSearch}
\alias{Jaatha.initialSearch}
\title{Search the parameter space for good starting positions}
\usage{
  Jaatha.initialSearch(jObject, nSim = 200,
    nBlocksPerPar = 3, multiple = F)
}
\arguments{
  \item{jObject}{The Jaatha settings (create with
  \code{\link{Jaatha.initialize}})}

  \item{nSim}{The number of simulations that are performed
  in each bin}

  \item{nBlocksPerPar}{The number of block per parameter.
  Will result in nPar^nBlocksPerPar blocks}

  \item{multiple}{If False only the best starting postions
  will be returned.  Otherwise it will be a list with all
  postitions.}
}
\value{
  The best starting positions or a list of all starting
  positions (see "multiple" above)
}
\description{
  This functions devides the parameter space in different
  parts (blocks). In each block, simulations for different
  parameter combinations are run to roughly predict the
  combination with the highest score (which is equivalent
  to the highest composite log likelihood). This points can
  later be used as starting positions for the secound
  estimation phase of Jaatha
  (\code{\link{Jaatha.refineSearch}})
}

