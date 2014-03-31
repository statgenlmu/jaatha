\docType{class}
\name{Jaatha-class}
\alias{Jaatha-class}
\title{The "Jaatha" S4 class saves the basic parameters for a Jaatha estimation
procedure}
\description{
Slots: \describe{ \item{simFunc}{Function used for
simulating} \item{par.ranges}{A nx2 matrix stating the
ranges for the n model parameters we want to estimate. The
first row gives the lower range for the parameters, the
second row the upper ranges. Each row stands for one
parameter and the row-names will be used as names for the
parameters.} \item{sum.stats}{The observed summary
statistics.} \item{seeds}{A set of random seeds. First one
to generate the other two.  The next one is for the initial
search, the last is for the refined search}
\item{cores}{The number of CPU cores to use for
simulations} \item{use.shm}{Use the shared memory /dev/shm
for temporary files (Linux only)} \item{opts}{Placeholder
for additional arguments, for instance ones that needs to
be passed to simFunc.} \item{calls}{The function calls for
the initial & refined search. These are used to rerun the
searches with the exactly same settings for generating
bootstrap confidence intervals and the likelihood-ratio
statistic.} \item{starting.positions}{A list of the
starting positions, returned by the initial search}
\item{likelihood.table}{A matrix with the best composite
log likelihood values and corresponding parameters}
\item{conf.ints}{Confidence Intervals for parameter
estimates produced by Jaatha.confidenceIntervals}
\item{route}{Tracks the best estimates of each step.} }
}

