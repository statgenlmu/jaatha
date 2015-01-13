#' @importFrom R6 R6Class
Stat_Base <- R6Class("Stat_Base",
  public = list(
    initialize = function(sim_data) self$set_data(sim_data),
    transform = function(sim_data) as.vector(sim_data$data),
    get_data = function() private$data,
    set_data = function(sim_data) private$data = self$transform(sim_data)
  ),
  private = list(
    data = NA
  )
)

#Stat_PoiInd <- R6Class("Stat_PoiInd", inherit = Stat_Base)

#' @importFrom reshape2 melt
Stat_PoiSmooth <- R6Class("Stat_PoiSmooth", inherit = Stat_Base,
  public = list(
    get_model = function() private$model,
    transform = function(sim_data) {
      data <- sim_data$data
      # Convert to data.frame
      dim_names <- lapply(dim(data), function(x) 1:x)
      names(dim_names) <- paste0('X', 1:length(dim(data)))
      dimnames(data) <- dim_names
      melt(data, value.name = 'sum.stat')
    },
    initialize = function(data, model) {
      private$data = self$transform(list(data=data))
      private$model = model
    }
  ),
  private = list(
    model = ""
  )
)