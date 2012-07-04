\name{dm.simSumStats}
\alias{dm.simSumStats}
\title{Simulates data according to a demographic model and calculates summary statistics form it}
\usage{
  dm.simSumStats(dm, parameters, sumStatFunc)
}
\arguments{
  \item{dm}{The demographic model according to which the
  simulations should be done}

  \item{parameters}{A vector of parameters which should be
  used for the simulations.  If a matrix is given, a
  simulation for each row of the matrix will be performed}

  \item{sumStatFunc}{The function to calculate summary
  statistics of the JSFS. It must return a numeric vector.}
}
\value{
  A matrix where each row is the vector of summary
  statistics for the parameters in the same row of the
  "parameter" matrix
}
\description{
  Simulates data according to a demographic model and
  calculates summary statistics form it
}
\examples{
dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100,seqLength=1000)
dm <- dm.addSpeciationEvent(dm,0.01,5)
dm <- dm.addMutation(dm,1,20)
dm.simSumStats(dm,c(1,10))
}

