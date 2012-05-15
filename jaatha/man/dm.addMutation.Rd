\name{dm.addMutation}
\alias{dm.addMutation}
\title{Adds mutations to a demographic model}
\usage{
  dm.addMutation(dm, lowerRange, upperRange, fixedValue)
}
\arguments{
  \item{dm}{The demographic model to which mutations should
  be added}

  \item{lowerRange}{If you want to estimate the mutation
  rate, this will be used as the smallest possible value.}

  \item{upperRange}{If you want to estimate the mutation
  rate, this will be used as the largest possible value.}

  \item{fixedValue}{If specified, the mutation rate will
  not be estimated, but assumend to be fixed at the given
  value.}
}
\value{
  The demographic model with mutation
}
\description{
  Adds mutations to a demographic model
}
\examples{
dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100,seqLength=1000)
dm <- dm.addSpeciationEvent(dm,0.01,5)
dm <- dm.addMutation(dm,1,20)
dm
}

