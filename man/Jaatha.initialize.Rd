\name{Jaatha.initialize}
\alias{Jaatha.initialize}
\title{Initialization of a Jaatha estimation for population genetics}
\usage{
  Jaatha.initialize(demographic.model, jsfs, seed,
    cores = 1, scaling.factor = 1, use.shm = FALSE,
    folded = FALSE, smoothing = FALSE)
}
\arguments{
  \item{demographic.model}{The demographic model to use}

  \item{jsfs}{Your observed Joint Site Frequency Spectrum
  (JSFS). Jaatha uses the JSFS as summary statistics.}

  \item{folded}{If 'TRUE', Jaatha will assume that the JSFS
  is folded.}

  \item{seed}{An integer used as seed for both Jaatha and
  the simulation software}

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
  Usually gives a huge performance boost.  Warning: This
  option will be removed in a future version of jaatha. The
  cleaner way to achieve this is to move your complete
  R-temp directory to the memory drive. This is explained
  on http://www.paulstaab.de/2013/11/r-shm .}

  \item{smoothing}{If set to true, Jaatha uses a different
  way to summaries the JSFS. Instead of binning certain
  areas, and fitting a glm per area, only one glm is fitted
  for the complete JSFS, and the position of the different
  entries is treated as a model parameter. This feature is
  still experimental and not recommended for productive use
  at the moment.}
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

