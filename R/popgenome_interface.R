#' Use PopGenome Data for Jaatha
#' 
#' This allows to use genetic data that was imported with the package
#' \pkg{PopGenome} to be analysed with Jaatha. As this relies on \pkg{coala}
#' to calculate the summary statistics, you need to use a model that was created
#' from a coala model.
#' 
#' @inheritParams create_jaatha_data
#' @param data The \code{GENOME} data from \pkg{PopGenome}.
#' @param coala_model The coala model that was used for creating \code{model}.
#' @param only_synonymous Only use synonymous SNP if set to \code{TRUE}. 
#'   This requires that \pkg{PopGenome} knows where coding regions are., e.g.
#'   by using gff files.
#' @param trios If you are using locus trios for your anaylsis, then you need
#'   to tell Jaatha with loci should be combined to a trio. To do so,
#'   set this parameter to a list, in which each entry is a numeric
#'   vector of length three containing the indexes of three loci 
#'   which will form a trio.
#' @export
create_jaatha_data.GENOME <- function(data, model, coala_model,
                                      only_synonymous = FALSE,
                                      trios = NULL,
                                      ...) {
  
  require_package("coala")
  
  check_popgenome_consistency(data, coala_model)
  seg_sites <- get_popgenome_segsites(data, only_synonymous, trios)
  
  sumstat <- lapply(coala::get_summary_statistics(coala_model), function(stat) {
    stat$calculate(seg_sites, NULL, NULL, coala_model)
  })
  
  create_jaatha_data(sumstat, model)
}


check_popgenome_consistency <- function(data, coala_model) {
  require_package("PopGenome")
  require_package("coala")
  
  # Check Populations
  for (pop in seq(along = data@populations)) {
    if (length(data@populations[[pop]]) != 
        length(coala::get_population_indiviuals(coala_model, pop))) {
      stop("Population ", pop, " has a different number of samples in the data",
           " and the coala model")
    }
    #if (coala::get_outgroup_size(coala_model) != length(data@outgroup)) {
    #  stop("Outgroup has a different number of samples in the data and the",
    #       " and the coala model")
    #}
  }
}
 
 
get_popgenome_segsites <- function(data, only_synonymous, trios) {
  require_package("PopGenome")

  if (is.null(trios)) {
    seg_sites_list <- lapply(1:length(data@n.valid.sites), function(i) {
      if (data@n.valid.sites[[i]] == 0) return(NULL)
      get_popgenome_locus(data, i, only_synonymous)
    })
  } else {
    assert_that(is.list(trios))
    seg_sites_list <- lapply(trios, function(trio) {
      assert_that(is.numeric(trio))
      assert_that(length(trio) == 3)
      left <- get_popgenome_locus(data, trio[1], only_synonymous)
      middle <- get_popgenome_locus(data, trio[2], only_synonymous)
      right <- get_popgenome_locus(data, trio[3], only_synonymous)
      
      seg_sites <- cbind(left, middle, right)
      
      attr(seg_sites, "positions") <- c(attr(left, "positions"),
                                        attr(middle, "positions"),
                                        attr(right, "positions"))
      
      attr(seg_sites, "locus") <- c(rep(-1, ncol(left)),
                                    rep( 0, ncol(middle)),
                                    rep( 1, ncol(right)))
      seg_sites
    })
  }
  
  seg_sites_list[!sapply(seg_sites_list, is.null)]
}


# Gets PopGenome's biallelic matrix (bam) and converts it to Jaatha's 
# segregating sites
get_popgenome_locus <- function(data, locus_number, only_synonymous) {
  require_package("PopGenome")
  
  bam <- PopGenome::get.biallelic.matrix(data, locus_number)
  if (is.null(bam)) {
    warning("Locus ", locus_number, " is NULL")
    seg_sites <- matrix(0, length(data@populations[[1]]) + 
                          length(data@populations[[2]]) + 
                          length(data@outgroup), 0)
    attr(seg_sites, "positions") <- numeric(0)
    return(seg_sites)
  }
  
  pop1 <- which(row.names(bam) %in% data@populations[[1]])
  pop2 <- which(row.names(bam) %in% data@populations[[2]])
  outgroup <- which(row.names(bam) %in% data@outgroup)
  stopifnot(!is.null(pop1))
  stopifnot(!is.null(pop2))
  stopifnot(!is.null(outgroup))
  
  # Select relevant data
  if (only_synonymous) {
    syn <- data@region.data@synonymous[[locus_number]]
    if (is.null(syn)) syn <- numeric(0)
    syn[is.na(syn)] <- FALSE
    seg_sites <- bam[c(pop1, pop2, outgroup), syn, drop = FALSE]
  } else {
    seg_sites <- bam[c(pop1, pop2, outgroup), , drop = FALSE]
  }
  stopifnot(!is.null(seg_sites))
  
  # Add positions attribute
  attr(seg_sites, "positions") <-
    as.numeric(colnames(bam)) / data@n.sites[[locus_number]]
  seg_sites
}


#' @importFrom utils capture.output
create_popgenome_test_data <- function() {
  require_package("PopGenome")
  
  # Create Test Data
  output <- tempfile("output")
  capture.output({
    fasta <- system.file("example_fasta_files", package = "jaatha")
    data_pg <- PopGenome::readData(fasta, progress_bar_switch = FALSE)
    data_pg <- PopGenome::set.outgroup(data_pg, c("Individual_Out-1", 
                                                  "Individual_Out-2"))
    data_pg <- PopGenome::set.populations(data_pg, 
                                          list(paste0("Individual_1-", 1:5), 
                                               paste0("Individual_2-", 1:5)))
  })
  data_pg
}
