\name{dm.addSymmetricMigration}
\alias{dm.addSymmetricMigration}
\title{Adds symmetric migration between all populations}
\usage{
  dm.addSymmetricMigration(dm, lower.range, upper.range,
    fixed.value, par.new = T, new.par.name = "M",
    parameter, time.start = "0")
}
\arguments{
  \item{dm}{The demographic model to which migration events
  should be added.}

  \item{par.new}{If 'TRUE' a new parameter will be created
  using the arguments 'lower.range' and 'upper.range' or
  'fixed.value'. It will be named 'new.par.name'.  If
  'FALSE' the argument 'parameter' will be evaluated
  instead.}

  \item{lower.range}{If you want to estimate the migration
  parameter (see note}

  \item{upper.range}{Same as lower.range, but the largest
  possible value.}

  \item{fixed.value}{If specified, the migration rate will
  not be estimated, but assumed to have the given value.}

  \item{new.par.name}{The name for the new parameter.}

  \item{parameter}{Instead of creating a new parameter, you
  can also set the mutation rate to an expression based on
  existing parameters. For example setting this to "M" will
  use an parameter with name M that you have previously
  created. You can also use R expression here, so "2*M" or
  "5*M+2*tau" (if tau is another parameter) will also work
  (also this does not make much sense).}

  \item{time.start}{The time point at which the migration
  with this rate starts.}
}
\value{
  The demographic model with migration
}
\description{
  This adds migration between all subpopulation, all with
  the same rate. Please look at the documentation for
  \link{dm.addMigration} for detailed information about
  migration.
}
\examples{
dm <- dm.createThetaTauModel(c(25,25), 100)
dm <- dm.addSymmetricMigration(dm, 0.01, 5, time.start="0.5*tau")
}

