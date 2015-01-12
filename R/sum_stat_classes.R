#' @importFrom R6 R6Class
Stat_Base <- R6Class("Stat_Base",
  public = list(
    transform = function(sim_data) as.vector(sim_data),
    get_data = function() private$data
  ),
  private = list(
    data = NA
  )
)

Stat_PoiInd <- R6Class("Stat_PoiInd", inherit = Stat_Base)

Stat_PoiSmooth <- R6Class("Stat_PoiSmooth", inherit = Stat_Base)