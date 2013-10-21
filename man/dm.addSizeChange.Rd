\name{dm.addSizeChange}
\alias{dm.addSizeChange}
\title{Adds an instantaneous change of the population size of one
population to a model.}
\usage{
  dm.addSizeChange(dm, min.size.factor, max.size.factor,
    fixed.size.factor, par.new = T, new.par.name = "q",
    parameter, population, at.time = "0")
}
\arguments{
  \item{dm}{The demographic model to which the size change
  should be added.}

  \item{par.new}{If 'TRUE' a new parameter will be created
  using the arguments 'min.size.factor' and
  'max.size.factor' or 'fixed.size.factor'. It will be
  named 'new.time.point.name' If 'FALSE' the argument
  'parameter' will be evaluated instead.}

  \item{min.size.factor}{If you want to estimate the size
  factor, this will be used as the smallest possible
  value.}

  \item{max.size.factor}{Same as min.size.factor, but the
  largest possible value.}

  \item{fixed.size.factor}{If specified, the size factor
  not be estimated, but assumed to have the given value.}

  \item{new.par.name}{Name for the new parameter.}

  \item{population}{The number of the population in which
  the spilt occurs. See \link{dm.addSpeciationEvent} for
  more information.}

  \item{parameter}{Instead of creating a new parameter, you
  can also set the mutation rate to an expression based on
  existing parameters. For example setting this to "tau"
  will use an parameter with name tau that you have
  previously created. You can also use R expression here,
  i.e. "2*tau" or "5*M+2*tau" (if M is another parameter)
  will also work (also this does not make much sense).}

  \item{at.time}{The time point at which the size changes.}
}
\value{
  The demographic model with a size change.
}
\description{
  This function changes the effective population size of
  one population. The change is performed at a given time
  point ('at.time') and applies to the time interval
  farther into the past from this point. The population
  size is set to a factor of the size of the ancestral
  population Ne.
}
\details{
  If you want to add a slow, continuous change over some
  time, then use the \link{dm.addGrowth} function.
}
\examples{
# A model with one smaller population
dm <- dm.createThetaTauModel(c(20,37), 88)
dm <- dm.addSizeChange(dm, 0.1, 1, population=2, at.time="0")
}

