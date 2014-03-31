\name{calculateJsfs}
\alias{calculateJsfs}
\title{Calculates the JSFS for data imported with ape}
\usage{
calculateJsfs(ape.data, pop1.rows, pop2.rows, outgroup.rows = NA)
}
\arguments{
  \item{ape.data}{DNA sequence data of multiple individuals
  from two populations and an optional outgroup sequence.
  It should be in a format as returned by ape's 'read.dna'
  function. Please refer to ape's manual for further
  details.}

  \item{pop1.rows}{A numeric vector indicating which
  individuals of your dataset belong to the first
  population. The individuals are revered by their position
  in the dataset/their row number in 'ape.data'.}

  \item{pop2.rows}{Same as 'pop1.rows', but for the
  individuals of the second population.}

  \item{outgroup.rows}{Same as 'pop1.rows', but for the
  individuals of the outgroup (if any). The outgroup can
  consist of more than one individual to account for
  ancestral misidentification. In this case, only positions
  in which all outgroup sequences are identically are
  considered.  If no outgroup is given a folded JSFS is
  calculated.}
}
\value{
The calculated JSFS, as matrix.
}
\description{
The function 'read.dna' of package 'ape' allows you do
import DNA data from various formats into R. This function
calculates the JSFS of such imported data which can then be
used as input for a Jaatha search. Please note that the
data must be aligned.
}

