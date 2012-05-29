\name{dm.addSpeciationEvent}
\alias{dm.addSpeciationEvent}
\title{Adds a speciation event to a demographic model}
\usage{
  dm.addSpeciationEvent(dm, lowerRange, upperRange,
    fixedValue)
}
\arguments{
  \item{dm}{The demographic model to which the speciation
  event should be added}

  \item{lowerRange}{If you want to estimate the mutation
  rate, this will be used as the smallest possible value.}

  \item{upperRange}{If you want to estimate the mutation
  rate, this will be used as the largest possible value.}

  \item{fixedValue}{If specified, the mutation rate will
  not be estimated, but assumend to be fixed at the given
  value.}
}
\value{
  The extended demographic model
}
\description{
  Adds a speciation event to a demographic model
}
\examples{
dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100,seqLength=1000)
dm <- dm.addSpeciationEvent(dm,0.01,5)
dm <- dm.addMutation(dm,1,20)
dm
}

