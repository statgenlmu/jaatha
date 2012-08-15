#!/usr/bin/Rscript --vanilla
#
# readFasta
# %DESCRIPTION%
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-08-14
# Licence:  GPLv3 or later
#

args <- commandArgs(TRUE)

# Enable just-in-time-compiler if availible
if ("compiler" %in% rownames(installed.packages())){
  library("compiler")
  invisible(compiler::enableJIT(3))
}

fasta2Jsfs <- function(population1, population2, outgroup) {
  pop1 <- readFasta(population1)
  pop2 <- readFasta(population2)
  outg <- readFasta(outgroup)
  
  mat <- cbind(outg, pop1, pop2)
  outg <- 1:ncol(outg)
  pop1 <- 1:ncol(pop1)
  pop2 <- 1:ncol(pop2)

  snp.mask <- apply(mat, 1, isSnp)
  snps  <- mat[snp.mask,]

  snp.types <- getSNPTypes(snps)
  jsfs <- calcJSFS(snp.types, c(length(pop1),length(pop2)))
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

isSnp <- function(gene.row) {
  # Filter out bad positions
  if (!all(gene.row %in% c("A","C","T","G"))) return(F)
  # Filter out disagreeing outgroups
  if (length(unique(gene.row[1:length(outg)])) != 1) return(F)
  # Filter out fixed positions
  if (all(gene.row == "A")) return(F)
  if (all(gene.row == "C")) return(F)
  if (all(gene.row == "T")) return(F)
  if (all(gene.row == "G")) return(F)
  return(T)
}

getSNPTypes <- function(snps){
  isSpanish <- c(rep(F, length(outg)),
                 rep(T, length(pop1)),
                 rep(F, length(pop2)))
  isItalian <- c(rep(F, length(outg)),
                 rep(F, length(pop1)),
                 rep(T, length(pop2)))

  snp.type <- matrix(0, nrow(snps), 2)

  for (i in 1:nrow(snps)){
    snp.row <- unlist(snps[i, ])
    derived <- snp.row != snp.row[1]
    snp.type[i, 1] <- sum(derived & isSpanish)
    snp.type[i, 2] <- sum(derived & isItalian)
  }
  return(snp.type)	
}

calcJSFS <- function(snps.list, sample.sizes) {
  jsfs <- matrix(0, sample.sizes[1] + 1, sample.sizes[2] + 1)

    for (j in 1:nrow(snp.types)) {
      jsfs <- jsfs + (row(jsfs) == (snp.types[j, 1] + 1) 
                      & col(jsfs) == (snp.types[j, 2] + 1))
    }
  return(jsfs)
}
