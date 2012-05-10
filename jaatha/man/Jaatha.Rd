\name{Jaatha}
\alias{Jaatha}
\title{The "Jaatha" S4 class saves the basic parameters for a Jaatha estimation
procedure}
\description{
  Slots: \describe{ \item{dm}{The demographic model to use}
  \item{nPar}{The number of parameters to estimate}
  \item{logMLmax}{The maximum log composite likelihood over
  all block} \item{nTotalSumstat}{The total number of
  summary statistics} \item{nNonJsfsSumstat}{The number of
  summary statistics that are not calculated from the JSFS}
  \item{popSampleSizes=}{Sample sizes of both populations.
  Vector of length 2.} \item{nLoci}{The number of loci}
  \item{MLest}{ML estimations of parameters (transformed
  between 0 and 1)} \item{seed}{Random seed for R}
  \item{externalTheta}{If TRUE theta will estimated using
  Waterson-estimator} \item{finiteSites}{If TRUE, we use a
  finite sites mutation model instead of an infinite sites
  one} \item{parNames}{The name of the paramters we
  estimate} \item{resultFile}{The name of a file in which
  the results will be written} \item{debugMode}{If TRUE, a
  debug output will be produced} \item{logFile}{If set, the
  debug output will be written into this file}
  \item{sumStats}{The summary statistics from the real
  data} \item{likelihood.table}{A matrix with the best
  composite log likelihood values and coresponding
  parameters} }
}

