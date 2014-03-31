\name{dm.addParameter}
\alias{dm.addParameter}
\title{Create a parameter that can be used for one or more features}
\usage{
dm.addParameter(dm, par.name, lower.boundary, upper.boundary, fixed.value)
}
\arguments{
  \item{dm}{The demographic model to which the parameter
  will be added}

  \item{par.name}{The name of the parameter. You can use
  this name later to access the parameter from a model
  feature}

  \item{lower.boundary}{The lower boundary of the range
  within which the parameter value will be estimated. Don't
  specify 'fixed.value' if you want to do so.}

  \item{upper.boundary}{Like 'lower.boundary', but the
  upper end of the parameter range.}

  \item{fixed.value}{If this argument is given, than rather
  than being estimated a fixed value will be used.}
}
\value{
The original model extended with the new parameter.
}
\description{
The function creates a new model parameter. It can either
have a fixed value or you can enter a range of possible
values if you want to estimate it.
}
\examples{
dm <- dm.createThetaTauModel(c(15,23), 100)
dm <- dm.addParameter(dm, "mig", 0.1, 5)
dm <- dm.addMigration(dm, par.new=FALSE, parameter="mig", pop.from=1, pop.to=2)
dm <- dm.addMigration(dm, par.new=FALSE, parameter="2*mig", pop.from=2, pop.to=1)
}

