\name{dm.defaultSumStats}
\alias{dm.defaultSumStats}
\title{These are the default summary statistics for Jaatha}
\usage{
  dm.defaultSumStats(jsfs)
}
\arguments{
  \item{jsfs}{The joint site frequency spectrum of two
  populations}
}
\value{
  A vector with sums over different areas of the JSFS
}
\description{
  These are the default summary statistics for Jaatha
}
\examples{
jsfs <- matrix(rpois(26*26,5),26,26)
dm.defaultSumStats(jsfs)
}

