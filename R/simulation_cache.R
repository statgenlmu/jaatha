#' @importFrom R6 R6Class
sim_cache_class <- R6Class("sim_cache",
  private = list(sim_data = list()),
  public = list(
    add = function(sim_data) {
      assert_that(is.list(sim_data))
      is_valid <- vapply(sim_data, function(x) {
        if (!is.list(x)) return(FALSE)
        if (!any(names(x) == "pars_normal")) return(FALSE)
        TRUE
      }, logical(1))
      if (!all(is_valid)) stop("Invalid simulation data")
      
      private$sim_data[1:length(sim_data) + length(private$sim_data)] <- 
        sim_data[1:length(sim_data)]
    },
    get_sim_data = function(block) {
      in_block <- vapply(private$sim_data, function(x) {
        block$includes(x$pars_normal)
      }, logical(1))
      private$sim_data[in_block]
    },
    get_size = function() length(private$sim_data)
  )
)

create_sim_cache <- function() {
  sim_cache_class$new()
}
