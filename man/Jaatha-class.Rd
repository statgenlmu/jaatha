\name{Jaatha-class}
\alias{Jaatha-class}
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
  between 0 and 1)} \item{seeds}{A set of random seeds.
  First one to generate the other two.  The next one is for
  the initial search, the last is for the refined search}
  \item{externalTheta}{Option removed in Version 2.1. Use a
  Version 2.0.2 if you want to use external Theta.}
  \item{finiteSites}{If TRUE, we use a finite sites
  mutation model instead of an infinite sites one}
  \item{parNames}{The name of the parameters we estimate}
  \item{debugMode}{If TRUE, a debug output will be
  produced} \item{logFile}{If set, the debug output will be
  written into this file} \item{sumStats}{The summary
  statistics from the real data} \item{likelihood.table}{A
  matrix with the best composite log likelihood values and
  corresponding parameters} \item{sum.stats.func}{The
  function for summarizing the JSFS}
  \item{sim.package.size}{Execute this number of
  simulations on a core in a row} \item{cores}{The number
  of CPU cores to use for simulations}
  \item{scaling.factor}{Only simulate this part of the data
  and interpolate the rest} \item{route}{Tracks the best
  estimates of each step.} \item{use.shm}{Use the shared
  memory /dev/shm for temporary files (Linux only)} }
}

