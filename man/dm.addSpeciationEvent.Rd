\name{dm.addSpeciationEvent}
\alias{dm.addSpeciationEvent}
\title{Adds a speciation event to a demographic model}
\usage{
  dm.addSpeciationEvent(dm, min.time, max.time, fixed.time,
    in.population = 1, new.time.point = T,
    new.time.point.name = NA, time.point = NA)
}
\arguments{
  \item{dm}{The demographic model to which the split should
  be added.}

  \item{new.time.point}{If 'TRUE' a new parameter will be
  created using the arguments 'min.time' and 'max.time' or
  'fixed.time'. It will be named 'new.time.point.name'.  If
  'FALSE' the argument 'time.point' will be evaluated
  instead.}

  \item{min.time}{If you want to estimate the time point,
  this will be used as the smallest possible value.}

  \item{max.time}{Same as min.time, but the largest
  possible value.}

  \item{fixed.time}{If specified, the time.point will not
  be estimated, but assumed to have the given value.}

  \item{new.time.point.name}{The name for the new time
  point.}

  \item{in.population}{The number of the population in
  which the spilt occurs. See above for more information.}

  \item{time.point}{Instead of creating a new parameter,
  you can also set the mutation rate to an expression based
  on existing parameters. For example setting this to "tau"
  will use an parameter with name tau that you have
  previously created. You can also use R expression here,
  i.e. "2*tau" or "5*M+2*tau" (if M is another parameter)
  will also work (also this does not make much sense).}
}
\value{
  The demographic model with a split.
}
\description{
  You can use this function the create a new population
  that splits of from an existing population at a given
  time in the past. The time can be given as parameter or
  as an expression based on previously generated time
  points.
}
\details{
  As always, time in measured in Number of 4Ne generations
  in the past, where Ne is the (effective) size of the
  ancestral population.

  The command will print the number of the new population,
  which will be the number of previously existing
  populations plus one. The ancestral population has number
  "1".
}
\examples{
dm <- dm.createDemographicModel(c(25,25), 100)
dm <- dm.addSpeciationEvent(dm,0.01,5)
dm <- dm.addMutation(dm,1,20)
}

