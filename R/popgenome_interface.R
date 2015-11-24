#' Use PopGenome Data for Jaatha
#' 
#' This allows to use genetic data that was imported with the package
#' \pkg{PopGenome} to be analysed with Jaatha. As this relies on \pkg{coala}
#' to calculate the summary statistics, you need to use a model that was created
#' from a coala model.
#' 
#' @inheritParams create_jaatha_data
#' @param data The \code{GENOME} data from \pkg{PopGenome}.

#' @param coala_model The coala model that was provided to 
#'   \code{\link{create_jaatha_model}}.
#' @param trios An optional argument that will be passed on to 
#'   \code{coala::calc_sumstats_from_data} and is documented there.
#' @param only_synonymous Only use synonymous SNP if set to \code{TRUE}. 
#'   This requires that \pkg{PopGenome} knows where coding regions are., e.g.
#'   by using gff files.
#' @export
create_jaatha_data.GENOME <- function(data, model, coala_model,
                                      trios = NULL,
                                      only_synonymous = FALSE,
                                      ...) {
  require_package("coala")
  check_popgenome_consistency(data, coala_model)
  
  segsites_list <- get_popgenome_segsites(data, only_synonymous)
  sumstats <- coala::calc_sumstats_from_data(coala_model, segsites_list)
  create_jaatha_data(sumstats, model, coala_model, trios, ...)
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
 
 
get_popgenome_segsites <- function(data, only_synonymous) {
  require_package("PopGenome")

  seg_sites_list <- lapply(seq_along(data@n.valid.sites), get_popgenome_locus,
                           data = data, only_synonymous = only_synonymous)
  
  seg_sites_list[!vapply(seg_sites_list, is.null, logical(1))]
}


# Gets PopGenome's biallelic matrix (bam) and converts it to Jaatha's 
# segregating sites
get_popgenome_locus <- function(locus_number, data, only_synonymous) {
  require_package("PopGenome")
  require_package("coala")
  
  if (data@n.valid.sites[[locus_number]] == 0) return(NULL)
  
  bam <- PopGenome::get.biallelic.matrix(data, locus_number)
  if (is.null(bam)) {
    warning("Locus ", locus_number, " is NULL")
    seg_sites <- coala::create_segsites(matrix(0, 
                                               length(data@populations[[1]]) + 
                                                 length(data@populations[[2]]) + 
                                                 length(data@outgroup), 0), 
                                        numeric(0))
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
  pos <- as.numeric(colnames(bam)) / data@n.sites[[locus_number]]
  coala::create_segsites(seg_sites, pos)
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
