\name{Jaatha.initialize}
\alias{Jaatha.initialize}
\title{Initialization of a Jaatha estimation for population genetics}
\usage{
  Jaatha.initialize(demographic.model, jsfs, seed,
    sim.package.size = 10, cores = 1, scaling.factor = 1,
    use.shm = FALSE, folded = FALSE)
}
\arguments{
  \item{demographic.model}{The demographic model to use}

  \item{jsfs}{Your observed Joint Site Frequency Spectrum
  (JSFS). Jaatha uses the JSFS as summary statistics.}

  \item{folded}{If 'TRUE', Jaatha will assume that the JSFS
  is folded.}

  \item{seed}{An integer used as seed for both Jaatha and
  the simulation software}

  \item{sim.package.size}{When running Jaatha on multiple
  cores, a singe core will always execute a whole "package"
  of simulations the reduce the inter thread communication
  overhead. This gives the number of simulations in such a
  package. Choose a number such that the execution of the
  simulations takes at least 15 seconds.}

  \item{cores}{The number of cores to use in parallel. If
  0, it tries to guess the number of available cores and
  use them all.}

  \item{scaling.factor}{You can use this option if you have
  a large dataset. If so, Jaatha only simulates only a
  fraction 1/scaling.factor of the dataset and interpolates
  the missing data.}

  \item{use.shm}{Logical. Many modern linux distributions
  have a shared memory file system available under
  /dev/shm. Set this to TRUE to use it for temporary files.
  Usually gives a huge performance boost.}
}
\value{
  A S4-Object of type jaatha containing the settings
}
\description{
  This function sets the basic parameters for an analysis
  with Jaatha and is the first step for each application of
  it.
}
\examples{
dm <- dm.createThetaTauModel(c(20,25), 100)
jsfs <- matrix(rpois(21*26, 5), 21, 26)
jaatha <- Jaatha.initialize(dm, jsfs)
}

