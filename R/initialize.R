#' Initialization of a Jaatha estimation for population genetics
#'
#' This function sets the basic parameters for an analysis with
#' Jaatha and is the first step for each application of it.
#'
#' @param model The demographic model to use
#' @param data  The observed data. Jaatha can use data imported with package 
#'              \pkg{PopGenome}. Please refer the to the vignette 
#'              "The Jaatha HowTo" for more information.
#' @param cores The number of cores to use in parallel. If 0, it tries to
#'              guess the number of available cores and use them all.
#' @param scaling_factor You can use this option if you have a large dataset. If
#'              so, Jaatha only simulates only a fraction 1/scaling_factor of the
#'              dataset and interpolates the missing data.
#' @param smoothing If set to true, Jaatha uses a different way to summaries the
#'              JSFS. Instead of binning certain areas, and fitting a glm per
#'              area, only one glm is fitted for the complete JSFS, and the
#'              position of the different entries is treated as a model
#'              parameter. This feature is still experimental and not
#'              recommended for productive use at the moment.
#' @param only_synonymous Only use synonymous SNP if set to \code{TRUE}. Requires
#'              to provided \code{data} as a PopGenome "GENOME" object.
#' @param trios If you are using locus trios for your anaylsis, then you need
#'              to tell Jaatha with loci should be combined to a trio. To do so,
#'              set this parameter to a list, in which each entry is a numeric
#'              vector of length three containing the indexes of three loci 
#'              which will form a trio.
#' @return A S4-Object of type jaatha containing the settings
#' @importFrom coala get_parameter_table get_summary_statistics 
#' @importFrom coala get_locus_number scale_model
#' @importFrom methods new representation
#' @importFrom assertthat assert_that
#' @export
Jaatha.initialize <- function(data, model, cores=1, scaling_factor=1,
                              smoothing=FALSE, only_synonymous=FALSE,
                              trios = NULL) {
  
  # --- Check parameters -------------------------------------
  assert_that("Coalmodel" %in% class(model)) 
  assert_that(is.numeric(cores))
  assert_that(length(cores) == 1)
  assert_that(is.numeric(scaling_factor))
  assert_that(length(scaling_factor) == 1)
  assert_that(is.logical(smoothing))
  assert_that(length(smoothing) == 1)
  assert_that(is.logical(only_synonymous))
  assert_that(length(only_synonymous) == 1)  
  
  # --- Convert the data into a list containing the seg.sites of the different groups
  if ("GENOME" %in% is(data)) {
    #checkModelDataConsistency(data, model)
    data <- convPopGenomeToSegSites(data, only_synonymous, trios)
  }
  if (!is.list(data)) stop("`data` has an unexpected format.")
  
  # ------------------------------------------------------------
  # Create Summary Statistics for summary statistic of the model
  # ------------------------------------------------------------
  sumstats <- list()
  
  model_sumstats <- get_summary_statistics(model)
  seg_sites <- data[["seg_sites"]]
  group <- 0
  if (is.null(seg_sites)) stop("No seg_sites in `data` for group ", group)
  assert_that(is.list(seg_sites))
  assert_that(length(seg_sites) == get_locus_number(model))
  assert_that(all(sapply(seg_sites, is.matrix)))
  
  for (sumstat in model_sumstats) {
    name <- sumstat$get_name()
    
    # --- JSFS Summary Statistic ------------------------------------
    if ("SumstatJsfs" %in% class(sumstat)) {
      if (!smoothing) {
        sumstats[[name]] <- Stat_JSFS$new(seg_sites, model, sumstat)
      } else {
        sumstats[[name]] <- 
          Stat_JSFS_smooth$new(seg_sites, model, sumstat)
        sumstats[[paste0("border_", name)]] <- 
          Stat_JSFS_border$new(seg_sites, model, sumstat)
      }
    }
    
    # --- JSFS Summary Statistic ------------------------------------
    else if ("SumstatSfs" %in% class(sumstat)) {
      sumstats[[name]] <- Stat_sfs$new(seg_sites, model, sumstat)
    }
    
    # --- Four Gamete Summary Statistic -----------------------------
    else if ("SumstatFourGamete" %in% class(sumstat)) {
      sumstats[[name]] <- Stat_FPC$new(seg_sites, model, sumstat)
    }
    
    # --- iHH Summary Statistic -------------------------------------
    else if ("sumstat_ihh" %in% class(sumstat)) {
      sumstats[[name]] <- Stat_Ihh$new(seg_sites, model, 
                                       sumstat, c(.25, .5, .75, .95))
    }
    
    # --- Omega" Summary Statistic ----------------------------------
    else if ("SumstatOmegaPrime" %in% class(sumstat)) {
      sumstats[[name]] <- Stat_OmegaPrime$new(seg_sites, model, 
                                              sumstat, c(.5, .75, .95))
    }
    
  }
  
  if (scaling_factor != 1) {
    model <- scale_model(model, scaling_factor)
  }
  
  
  # ------------------------------------------------------------
  # Create the Jaatha object
  # ------------------------------------------------------------
  par_ranges <- as.matrix(get_parameter_table(model)[,-1])
  rownames(par_ranges) <- get_parameter_table(model)$name
  
  jaatha <- new("Jaatha", 
                sim_func=function(sim.pars, jaatha) {
                  simulate(jaatha@opts[["model"]], pars=sim.pars)
                },
                par_ranges=par_ranges,  
                sum_stats=sumstats,
                cores=cores,
                options = list(model=model),
                scaling_factor = scaling_factor)
  
  invisible(jaatha)
}