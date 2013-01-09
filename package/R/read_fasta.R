#-------------------------------------------------------------------------------
# read_fasta.R
# Methodes for calculating the JSFS of two species based on (aligned) sequences
# stored in fasta format.
# 
# Author:   Paul R. Staab & Lisha Naduvilezhath 
# Email:    staab (at) bio.lmu.de
# Date:     2012-09-20
# Licence:  GPLv3 or later
#-------------------------------------------------------------------------------

#' Calculates the JSFS for data imported with ape
#'
#' The function 'read.dna' of package 'ape' allows you do import DNA data from
#' various formats into R. This function calculates the JSFS of such imported
#' data which can then be used as input for a Jaatha search. Please note that
#' the data must be aligned.
#'
#' @param ape.data    DNA sequence data of multiple individuals from two populations
#'                    and an optional outgroup sequence. It should be in a
#'                    format as returned by ape's 'read.dna' function. Please
#'                    refer to ape's manual for further details.
#' @param pop1.rows   A numeric vector indicating which individuals of your
#'                    dataset belong to the first population. The individuals
#'                    are revered by their position in the dataset/their row
#'                    number in 'ape.data'.
#' @param pop2.rows   Same as 'pop1.rows', but for the individuals of the second
#'                    population. 
#' @param outgroup.rows  Same as 'pop1.rows', but for the individuals of the
#'                    outgroup (if any). The outgroup can consist of more than
#'                    one individual to account for ancestral
#'                    misidentification. In this case, only positions in which
#'                    all outgroup sequences are identically are considered.
#'                    If no outgroup is given a folded JSFS is calculated.
#' @return            The calculated JSFS, as matrix.
#' @export
calculateJsfs <- function(ape.data, pop1.rows, pop2.rows, outgroup.rows=NA) {
  jaatha:::checkType(pop1.rows, c("num", "vec"), T, F)
  jaatha:::checkType(pop2.rows, c("num", "vec"), T, F)
  jaatha:::checkType(outgroup.rows, c("num", "vec"), F, T)

  if (all(is(ape.data) == "DNAbin")) ape.data <- as.character(ape.data)
  jaatha:::checkType(ape.data, c("mat"), F, T)

  jsfs <- markerTableToJsfs(t(ape.data), pop1.rows, pop2.rows, outgroup.rows)
  return(jsfs)
}


markerTableToJsfs <- function(marker.table, pop1.cols, pop2.cols,
                              outgroup.cols=NA){
  
  folded <- F
  if (any(is.na(outgroup.cols))) folded <- T

  pop1.mask <- 1:ncol(marker.table) %in% pop1.cols
  pop2.mask <- 1:ncol(marker.table) %in% pop2.cols
  outg.mask <- 1:ncol(marker.table) %in% outgroup.cols
 
  snp.mask <- apply(marker.table, 1, isSnp, outg.mask=outg.mask)
  snps  <- marker.table[snp.mask,]

  snp.types <- getSNPTypes(snps, pop1.mask, pop2.mask, folded)
  jsfs <- calcJSFS(snp.types, c(length(pop1.cols),length(pop2.cols)))
  return(jsfs)
}


isSnp <- function(gene.row, outg.mask) {
  # Filter out bad positions
  if (!all(gene.row %in% c("a","c","t","g"))) return(F)
  # Filter out fixed positions
  if (all(gene.row == "A")) return(F)
  if (all(gene.row == "C")) return(F)
  if (all(gene.row == "T")) return(F)
  if (all(gene.row == "G")) return(F)
  # Now we are finished if we have no outgrout
  if (all(!outg.mask)) {
    if (length(unique(gene.row)) > 2) return(F)
    return(T)
  }
  # Otherwise filter out disagreeing outgroups
  if (length(unique(gene.row[outg.mask])) != 1) return(F)
  return(T)
}


getSNPTypes <- function(snps, pop1.mask, pop2.mask, folded){
  snp.type <- matrix(0, nrow(snps), 2)
  pops.size <- c(sum(pop1.mask), sum(pop2.mask))

  if (nrow(snps) == 0) return(snp.type)

  for (i in 1:nrow(snps)){
    snp.row <- unlist(snps[i, ])
    derived <- snp.row != snp.row[1]
    snp.type[i, 1] <- sum(derived & pop1.mask)
    snp.type[i, 2] <- sum(derived & pop2.mask)
    if (folded) {
      if (sum(snp.type[i, ]) > ceiling(sum(pops.size)/2))
        snp.type[i, ] <- pops.size - snp.type[i, ]
    }
  }
  return(snp.type)	
}


calcJSFS <- function(snp.types, sample.sizes) {
  jsfs <- matrix(0, sample.sizes[1] + 1, sample.sizes[2] + 1)
  
  if (nrow(snp.types) == 0) return(jsfs)

    for (j in 1:nrow(snp.types)) {
      jsfs <- jsfs + (row(jsfs) == (snp.types[j, 1] + 1) 
                      & col(jsfs) == (snp.types[j, 2] + 1))
    }
  return(jsfs)
}
