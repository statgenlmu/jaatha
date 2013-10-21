\name{addFeature}
\alias{addFeature}
\title{Low level function for adding a new feature to a demographic Model}
\usage{
  addFeature(dm, type, parameter = NA, lower.range = NA,
    upper.range = NA, fixed.value = NA, par.new = T,
    pop.source = NA, pop.sink = NA, time.point = NA,
    group = NA)
}
\arguments{
  \item{dm}{The demographic model to which the feature will
  be added}

  \item{type}{The type of the feature coded as a character}

  \item{lower.range}{The lower boundary for the value of
  the parameter}

  \item{upper.range}{The upper boundary for the value of
  the parameter}

  \item{fixed.value}{If given, the parameter will be set to
  a fixed value. This is equivalent to seting lower.range
  equal to upper.range.}

  \item{pop.source}{The source population if availible
  (think e.g. of migration)}

  \item{pop.sink}{The target or "sink" population if
  availible (think e.g. of migration)}

  \item{time.point}{Normally the point in backwards time
  where the feature starts.}

  \item{par.new}{If 'TRUE' a new parameter will be created
  using the arguments 'lower.range' and 'upper.range' or
  'fixed.value'. It will be named 'parameter'.  If 'FALSE'
  the argument 'parameter' will be evaluated instead.}

  \item{parameter}{Either the name of the parameter
  (par.new=TRUE), or an R expression possibly containing
  one or more previously created parameter names.}

  \item{group}{For genomic features, different groups can
  be created.}
}
\value{
  The extended demographic model.
}
\description{
  Use this function to add a feature to a dm if there is no
  higher level "dm.add*"-function availible.
}

