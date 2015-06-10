#' @importFrom R6 R6Class
Stat_Cube <- R6Class("Stat_Cube", inherit = Stat_PoiInd,
  private = list(
    stat_name = NA,
    breaks = NA
  ),
  public = list(
    to_matrix = function(x) x,
    initialize = function(seg_sites, model, stat, break_probs = c(.1, .5)) {
      if (any(break_probs > 1)) stop("probs greater then one")
                        
      private$stat_name <- stat$get_name()
      value <- stat$calculate(seg_sites, NULL, model)
      fake_sim_data <- list()
      fake_sim_data[[private$stat_name]] <- value
      
      value_matrix <- self$to_matrix(value)
      private$breaks = lapply(1:ncol(value_matrix), function(i) {
        calcBreaks(value_matrix[, i], break_probs)
      })
      names(private$breaks) <- colnames(value_matrix)
                        
      super$initialize(fake_sim_data, stat$get_name())
    },
    transform = function(sim_data) {
      stat <- self$to_matrix(sim_data[[private$stat_name]])                  
      generateLociCube(stat, private$breaks, 1:ncol(stat))
    },
    get_breaks = function() private$breaks
   )
)

#' @importFrom R6 R6Class
Stat_Ihh <- R6Class("Stat_Ihh", inherit = Stat_Cube,
  public = list(
    to_matrix = function(value) {
      do.call(rbind, lapply(value, function(x) max(x[ , 3])))
    }
  )
)


#' @importFrom R6 R6Class
Stat_OmegaPrime <- R6Class("Stat_OmegaPrime", inherit = Stat_Cube,
  public = list(
    to_matrix = function(value) {
      matrix(value, ncol = 1)
    }
  )
)
    
