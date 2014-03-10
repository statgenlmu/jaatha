\name{dm.addSampleSize}
\alias{dm.addSampleSize}
\title{Sets how many individuals from each population are sampled at time 0.}
\usage{
dm.addSampleSize(dm, sample.size, group = 0)
}
\arguments{
  \item{dm}{The demographic model to which recombination
  events should be added.}

  \item{sample.size}{A vector with sample sizes for each
  population.}

  \item{group}{The group of loci with this sample size.}
}
\value{
The demographic model with the sample
}
\description{
Sets how many individuals from each population are sampled
at time 0.
}

