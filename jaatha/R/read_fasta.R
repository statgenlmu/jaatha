#-------------------------------------------------------------------------------
# read_fasta.R
# Methodes for calculating the JSFS of two species based on (aligned) sequences
# stored in fasta format.
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-09-20
# Licence:  GPLv3 or later
#-------------------------------------------------------------------------------

#' Calculates the JSFS out of fasta files
#'
#' This function parses the three fasta file provided as input - one for each 
#' population and one for the outgroup - and calculates the Joint Site 
#' Frequency Spectrum (JSFS) of the sequences.
#' @param population1 The path to the fasta file which contains the aligned 
#'                    sequences of the first population.
#' @param population2 The path to the fasta file which contains the aligned 
#'                    sequences of the second population.
#' @param population2 The path to the fasta file which contains the aligned 
#'                    sequences of outgroup sequences. This file can contain 
#'                    multiple sequences to account for ancestral
#'                    misidentification. In this case, only positions in which
#'                    all outgroup sequences are identically are considered.
#' @return            The JSFS, as matrix.
#' @export
fasta2Jsfs <- function(population1, population2, outgroup) {
  pop1 <- readFasta(population1)
  pop2 <- readFasta(population2)
  outg <- readFasta(outgroup)
  
  mat <- cbind(outg, pop1, pop2)
  outg.cols <- 1:ncol(outg)
  pop1.cols <- 1:ncol(pop1) + ncol(outg)
  pop2.cols <- 1:ncol(pop2) + ncol(outg) + ncol(pop1)

  jsfs <- markerTableToJsfs(mat, pop1.cols, pop2.cols, outg.cols)
  return(jsfs)
}


markerTableToJsfs <- function(marker.table, pop1.cols, pop2.cols,
                              outgroup.cols){

  pop1.mask <- 1:ncol(marker.table) %in% pop1.cols
  pop2.mask <- 1:ncol(marker.table) %in% pop2.cols
  outg.mask <- 1:ncol(marker.table) %in% outgroup.cols

  snp.mask <- apply(marker.table, 1, isSnp, outg.mask=outg.mask)
  snps  <- marker.table[snp.mask,]

  snp.types <- getSNPTypes(snps, pop1.mask, pop2.mask)
  jsfs <- calcJSFS(snp.types, c(length(pop1.cols),length(pop2.cols)))
  return(jsfs)
}


readFasta <- function(fasta.file) {
  seq.data <- scan(fasta.file, character())
  first.char <- substr(seq.data,1,1)

  sequences <- list()
  individual <- 0

  for (i in seq(along = seq.data)) {
    if (first.char[i] == ">") {
      individual <- individual + 1
      sequences[[individual]] <- character()
    }
    else {
      seq.vector <- unlist(strsplit(seq.data[i], ""))
      sequences[[individual]] <- c(sequences[[individual]], seq.vector)
    }
  }
  return(sapply(sequences, function(x){return(x)}))
}


isSnp <- function(gene.row, outg.mask) {
  # Filter out bad positions
  if (!all(gene.row %in% c("A","C","T","G"))) return(F)
  # Filter out disagreeing outgroups
  if (length(unique(gene.row[outg.mask])) != 1) return(F)
  # Filter out fixed positions
  if (all(gene.row == "A")) return(F)
  if (all(gene.row == "C")) return(F)
  if (all(gene.row == "T")) return(F)
  if (all(gene.row == "G")) return(F)
  return(T)
}


getSNPTypes <- function(snps, pop1.mask, pop2.mask){
  snp.type <- matrix(0, nrow(snps), 2)

  for (i in 1:nrow(snps)){
    snp.row <- unlist(snps[i, ])
    derived <- snp.row != snp.row[1]
    snp.type[i, 1] <- sum(derived & pop1.mask)
    snp.type[i, 2] <- sum(derived & pop2.mask)
  }
  return(snp.type)	
}


calcJSFS <- function(snp.types, sample.sizes) {
  jsfs <- matrix(0, sample.sizes[1] + 1, sample.sizes[2] + 1)

    for (j in 1:nrow(snp.types)) {
      jsfs <- jsfs + (row(jsfs) == (snp.types[j, 1] + 1) 
                      & col(jsfs) == (snp.types[j, 2] + 1))
    }
  return(jsfs)
}
