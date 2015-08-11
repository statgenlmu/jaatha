jaatha_log_class <- R6Class("jaatha_log", 
  private = list(
    estimates = matrix(0, 0, 0),
    max_steps = 0,
    verbose = FALSE,
    par_ranges = NULL,
    format_par = function(par) {
      paste(format(private$par_ranges$denormalize(par)), collapse = " ")
    }
  ),
  public = list(
    initialize = function(reps, max_steps, model, verbose = TRUE) {
      par_number <- model$get_par_number()
      par_names <- model$get_par_ranges()$get_par_names()
      private$estimates <- matrix(NA, reps * max_steps, par_number + 3)
      colnames(private$estimates) <- c("rep", "step", "lh", par_names)
      private$max_steps <- max_steps
      private$verbose <- verbose
      private$par_ranges <- model$get_par_ranges()
    },
    log_estimate = function(rep, step, estimate) {
      idx <- (rep - 1) * private$max_steps + step
      assert_that(is_positive_int(idx))
      private$estimates[idx, ] <- c(rep, step, estimate$value, estimate$par)
      
      if (private$verbose) {
        message("Step ", step, ": Loglikelihood ", format(estimate$value), 
                ", Parameter: ", private$format_par(estimate$par))
      }
    },
    log_new_rep = function(rep, start_pos) {
      if (private$verbose) {
        message("Repetition ", rep, 
                " starting at ", private$format_par(start_pos))
      }
    },
    log_convergence = function() {
      if (private$verbose) message("Convergence detected")
    },
    get_estimates = function() private$estimates,
    create_results = function() {}
  )
)

create_jaatha_log <- jaatha_log_class$new
