\name{setLogging}
\alias{setLogging}
\title{Function to activate/change logging}
\usage{
  setLogging(log.level, log.file)
}
\arguments{
  \item{log.level}{An integer indicating the level of the
  log output.  If the level is smaller or equal to the log
  level stored in 'log.level', than the output is printed,
  otherwise it is discarded. 0 is the default level.}

  \item{log.file}{If given, the logging output will be
  written into this file.}
}
\value{
  nothing
}
\description{
  Function to activate/change logging
}

