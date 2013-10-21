\name{checkType}
\alias{checkType}
\title{Checks if a variable is of a given type and calls stop() on type mismatch}
\usage{
  checkType(variable, type, required = T, allow.na = T)
}
\arguments{
  \item{variable}{the variable to check}

  \item{type}{the name of the type the variable should
  have. Can be num/numeric, vec/vector mat/matrix or
  fun/function. Can also be s/single, in that case it must
  be a vector of length one .  If a vector of type names is
  given, the variable must be of all types.}

  \item{required}{A boolean that indicates whether the
  variable must be specified or can be missing. Value of
  NULL also counts as missing.}

  \item{allow.na}{If FALSE, an error is returned if the
  variable contains NA values. Only works for vectors.}
}
\value{
  nothing
}
\description{
  See heading. Side-effect warning: Calls stop() on type
  mismatch.
}

