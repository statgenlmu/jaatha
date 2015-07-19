#' @importFrom coala get_outgroup_size
checkModelDataConsistency <- function(data, model) {
  if (!requireNamespace("PopGenome", quietly = TRUE)) {
    stop("Please install package 'PopGenome'")
  }
  
  if (get_outgroup_size(model) != length(data@outgroup)) {
    stop("Expecting an outgroup size of ", length(data@outgroup), 
         " from data, but is ", get_outgroup_size(model), " in model.")
  }
}


convPopGenomeToSegSites <- function(data, 
                                    only_synonymous = FALSE, 
                                    trios = NULL) {
  
  if (is.null(trios)) {
    seg_sites_list <- lapply(1:length(data@n.valid.sites), function(i) {
      if (data@n.valid.sites[[i]] == 0) return(NULL)
      get_segsites(data, i, only_synonymous)
    })
  } else {
    assert_that(is.list(trios))
    seg_sites_list <- lapply(trios, function(trio) {
      assert_that(is.numeric(trio))
      assert_that(length(trio) == 3)
      left <- get_segsites(data, trio[1], only_synonymous)
      middle <- get_segsites(data, trio[2], only_synonymous)
      right <- get_segsites(data, trio[3], only_synonymous)
      
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

  list(seg_sites = seg_sites_list[!sapply(seg_sites_list, is.null)])
}


# Gets PopGenome's biallelic matrix (bam) and converts it to Jaatha's 
# segregating sites
get_segsites <- function(data, locus_number, only_synonymous) {
  if (!requireNamespace("PopGenome", quietly = TRUE)) {
    stop("Please install package 'PopGenome'")
  }
  
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


#' @importFrom coala coal_model feat_outgroup
createModelFromPopGenome <- function(data, quiet=FALSE) {
  if (!requireNamespace("PopGenome", quietly = TRUE)) {
    stop("Please install package 'PopGenome'")
  }
  
  stopifnot("GENOME" %in% is(data))
  sample_sizes <- sapply(data@populations, length)
  if (!quiet) { 
    message("Sample Sizes: ")
    for (pop in seq(along = sample_sizes)) {
      message(" Population ", pop, ": ", sample_sizes[pop])
    }
  }
  
  outgroup_size <- length(data@outgroup)
  outgroup_number <- length(sample_sizes) + 1
  if (!quiet) message("Outgroup size: ", paste(outgroup_size, collapse=" "), 
                      " (will be population ", outgroup_number, ")")


  loci_mask <- data@n.valid.sites > 0
  loci_length <- round(mean(data@n.valid.sites[loci_mask]))
  loci_number <- sum(loci_mask)
  
  if (!quiet) message("Number of Loci: ", loci_number)
  if (!quiet) message("Average Loci Length: ", loci_length, "bp")
  
  # Calculate TS/TV, but don"t add it to the model
  tstv_ratio <- sum(sapply(data@region.data@transitions[loci_mask], sum)) /
      sum(sapply(data@region.data@transitions[loci_mask], function(x) sum(1-x)))
  if (!quiet) message("Observed TS/TV: ", tstv_ratio, " (Not added to Model)")
  
  coal_model(c(sample_sizes, outgroup_size), loci_number, loci_length) +
    feat_outgroup(outgroup_number)
}
