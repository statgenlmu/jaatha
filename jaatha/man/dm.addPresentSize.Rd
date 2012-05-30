\name{dm.addPresentSize}
\alias{dm.addPresentSize}
\title{Adds a different present day population size of one population to a demographic model}
\usage{
  dm.addPresentSize(dm, lowerRange, upperRange, fixedValue,
    population)
}
\arguments{
  \item{dm}{The demographic model to which the speciation
  event should be added}

  \item{lowerRange}{If you want to estimate the size
  parameter, this will be used as the smallest possible
  value.}

  \item{upperRange}{If you want to estimate the size
  parameter, this will be used as the largest possible
  value.}

  \item{fixedValue}{If specified, the size parameter will
  not be estimated, but assumend to be fixed at the given
  value.}

  \item{population}{The population of which the size should
  change. Can be 1 or 2 at the moment.}
}
\value{
  The extended demographic model
}
\description{
  Adds a different present day population size of one
  population to a demographic model
}
\examples{
dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100)
dm <- dm.addSpeciationEvent(dm,0.01,5)
dm <- dm.addPresentSize(dm,2,3,population=2)
dm <- dm.addMutation(dm,1,20)
dm
}

