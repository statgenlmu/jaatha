\name{runSimulations}
\alias{runSimulations}
\title{The function that actually executes the simulations from within Jaatha.}
\usage{
runSimulations(pars, cores, jaatha)
}
\arguments{
  \item{pars}{A matrix with parameter values for the
  simulations. Each row is a set of model parameters.
  Should be in jaatha internal 0-1 scaling.}

  \item{cores}{The number of cores to use. Using more than
  one core requires the 'multicore' package.}

  \item{jaatha}{The current jaatha object}
}
\value{
A list, where each entry is a list of summary statistics
for a simulation.
}
\description{
The function that actually executes the simulations from
within Jaatha.
}

