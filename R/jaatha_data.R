#' @importFrom R6 R6Class
jaatha_data_class <- R6Class("jaatha_data", 
  lock_objects = TRUE, lock_class = TRUE,
  private = list(
    values = list(),
    options = list()
  ),
  public = list(
    initialize = function(data, model) {
      assert_that(is_jaatha_model(model))
      private$values <- lapply(model$get_sum_stats(), function(stat) {
        private$options[[stat$get_name()]] <- stat$generate_data_opts(data)
        stat$calculate(data, private$options[[stat$get_name()]])
      })
    },
    get_values = function(stat = NULL) {
      if (is.null(stat)) return(private$values)
      if (is.character(stat) || is.numeric(stat)) return(private$values[[stat]])
      private$values[[stat$get_name()]]
    },
    get_options = function(stat = NULL) {
      if (is.null(stat)) return(private$options)
      if (is.character(stat) || is.numeric(stat)) {
        return(private$options[[stat]])
      }
      private$options[[stat$get_name()]]
    }
  )
)


is_jaatha_data <- function(x) inherits(x, "jaatha_data")


create_jaatha_data <- function(data, model, ...) UseMethod("create_jaatha_data")


create_jaatha_data.default <- function(data, model, ...) {
  jaatha_data_class$new(data, model)
}


create_test_data <- function(model) {
  test_data <- model$test(quiet = TRUE)
  create_jaatha_data(test_data, model)
}
