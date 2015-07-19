#' @importFrom R6 R6Class
stat_basic_class <- R6Class("jaatha_stat_basic",
  lock_objects = FALSE, lock_class = TRUE,
  public = list(
    initialize = function(name, calc_func) {
      assert_that(is.character(name) && length(name) == 1)
      private$name <- name
      
      assert_that(is.function(calc_func))
      self$calculate <- calc_func
    },
    get_name = function() private$name,
    generate_data_opts = function(data) NULL
  ),
  private = list(
    name = ""
  )
)


create_jaatha_stat <- function(name, calc_func) {
  stat_basic_class$new(name, calc_func)
}


stat_identity <- function() create_jaatha_stat("id", function(x, y) x)
stat_sum <- function() create_jaatha_stat("sum", function(x, y) sum(x))


#' @importFrom reshape2 melt
Stat_PoiSmooth <- R6Class("Stat_PoiSmooth", inherit = stat_basic_class,
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
      names(dim_names) <- paste0("X", 1:length(dim(data)))
      dimnames(data) <- dim_names
      melt(data, value.name = "sum.stat")
    }
  )
)
