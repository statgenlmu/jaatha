#' @importFrom R6 R6Class
Stat_Base <- R6Class("Stat_Base",
  public = list(
    initialize = function(sim_data, name) {
      private$name = name
      self$set_data(sim_data)
      if (length(private$data) == 0) 
        stop('Failed to convert real data to ', name)
      if (any(is.na(private$data)))
        stop('NAs in real data for ', name)
    },
    transform = function(sim_data) as.vector(sim_data$data),
    get_data = function() private$data,
    set_data = function(sim_data) private$data = self$transform(sim_data),
    get_name = function() private$name
  ),
  private = list(
    data = NA,
    name = ''
  )
)

Stat_PoiInd <- R6Class("Stat_PoiInd", inherit = Stat_Base)

#' @importFrom reshape2 melt
Stat_PoiSmooth <- R6Class("Stat_PoiSmooth", inherit = Stat_Base,
  public = list(
    get_model = function() private$model,
    transform = function(sim_data) {
      private$to_data_frame(sim_data$data)
    },
    initialize = function(data, name, model) {
      super$initialize(data, name)
      private$model = model
    }
  ),
  private = list(
    model = "",
    to_data_frame = function(data) {
      stopifnot(!is.null(data))
      dim_names <- lapply(dim(data), function(x) 1:x)
      names(dim_names) <- paste0('X', 1:length(dim(data)))
      dimnames(data) <- dim_names
      melt(data, value.name = 'sum.stat')
    }
  )
)