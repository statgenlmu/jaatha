\name{dm.addRecombination}
\alias{dm.addRecombination}
\title{Adds recombination events to a demographic model}
\usage{
dm.addRecombination(dm, lower.range, upper.range, fixed.value, par.new = T,
  new.par.name = "rho", parameter)
}
\arguments{
  \item{dm}{The demographic model to which recombination
  events should be added.}

  \item{par.new}{If 'TRUE' a new parameter will be created
  using the arguments 'lower.range' and 'upper.range' or
  'fixed.value'. It will be named 'new.par.name'.  If
  'FALSE' the argument 'parameter' will be evaluated
  instead.}

  \item{lower.range}{If you want to estimate the
  recombination rate (see note}

  \item{upper.range}{Same as lower.range, but the largest
  possible value.}

  \item{fixed.value}{If specified, the mutation rate will
  not be estimated, but assumed to have the given value.}

  \item{new.par.name}{The name for the new parameter.}

  \item{parameter}{Instead of creating a new parameter, you
  can also set the mutation rate to an expression based on
  existing parameters. For example setting this to "rho"
  will use an parameter with name theta that you have
  previously created. You can also use R expression here,
  so "2*rho" or "5*rho+2*tau" (if tau is another parameter)
  will also work (also it does not make much sense).}
}
\value{
The demographic model with recombination
}
\description{
This function add the assumption to the model that
recombination events may occur within each locus. The
corresponding parameter - usually name rho - equals 4*Ne*r,
where Ne is the number of diploid individuals in the
ancestral population and r is the probability that a
recombination event within the locus will occur in one
generation. Even when using an infinite sites mutation
model, this assumes an finite locus length which is given
by the 'seq.length' parameter of the demographic model.
}
\details{
Please note that it does not make sense to estimate
recombination rates with Jaatha because it assumes unlinked
loci.
}
\examples{
dm <- dm.createDemographicModel(c(25,25), 100)
dm <- dm.addSpeciationEvent(dm, 0.01, 5)
dm <- dm.addRecombination(dm, fixed=20)
dm <- dm.addMutation(dm, 1, 20)
}

