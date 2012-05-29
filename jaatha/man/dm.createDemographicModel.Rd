\name{dm.createDemographicModel}
\alias{dm.createDemographicModel}
\title{Create a basic demographic model}
\usage{
  dm.createDemographicModel(sampleSizes, nLoci,
    seqLength = 1000, finiteSites = F, tsTvRatio = 0.33,
    debugMode = F, logFile = "")
}
\arguments{
  \item{sampleSizes}{Number of individual that are sampled.
  If your model has multiple populations, this needs to be
  a vector containing the sample sizes from each
  population.}

  \item{nLoci}{Number of loci that will be simulated}

  \item{seqLength}{Number of bases for each locus}

  \item{finiteSites}{If 'TRUE', a finite sites mutation
  model is assumed instead of an infinite sites one.}

  \item{tsTvRatio}{Transition transversion ratio}

  \item{debugMode}{If 'TRUE', a debug output will be
  produced}

  \item{logFile}{If set, the debug output will be written
  into the given file}
}
\value{
  The demographic model
}
\description{
  This function creates a basic empty demographic model,
  which is returned. Features like mutation, population
  splits and migration can be added afterwards.
}
\examples{
dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100,seqLength=1000)
dm <- dm.addSpeciationEvent(dm,0.01,5)
dm <- dm.addMutation(dm,1,20)
dm
}

