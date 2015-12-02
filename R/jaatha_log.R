jaatha_log_class <- R6Class("jaatha_log", 
  private = list(
    estimates = NULL,
    final_estimates = NULL,
    reps = 0,
    sim = 0,
    max_steps = 0,
    init_method = "none",
    verbose = FALSE,
    par_ranges = NULL,
    converged = NULL,
    format_par = function(par) {
      paste(format(private$par_ranges$denormalize(par)), collapse = " ")
    }
  ),
  public = list(
    initialize = function(model, data, reps, sim, max_steps, init_method, 
                          verbose = TRUE) {
      
      par_number <- model$get_par_number()
      par_names <- model$get_par_ranges()$get_par_names()
      private$estimates <- lapply(1:reps, function(i) {
        estimates <- matrix(NA, max_steps, par_number + 3)
        colnames(estimates) <- c("rep", "step", "llh", par_names)
        as.data.frame(estimates)
      })
      private$final_estimates <- private$estimates[[1]][rep(1, 5 * par_number), ]
      private$reps <- reps
      private$sim <- sim
      private$max_steps <- max_steps
      private$init_method <- init_method
      private$verbose <- verbose
      private$par_ranges <- model$get_par_ranges()
      private$converged <- rep(FALSE, reps)
    },
    log_estimate = function(rep, step, estimate, old_llh = NULL) {
      if (rep == "final") rep <- 0.0
      entry <- c(rep, step, estimate$value, estimate$par)
      if (rep == 0) {
        private$final_estimates[step, ] <- entry
      } else if (is.numeric(rep)) {
        private$estimates[[rep]][step, ] <- entry
      } else {
        stop("Unexpected value for 'rep'")
      }
      
      if (!private$verbose) return(invisible(NULL))
      if (rep != 0.0) {
        message(" Step ", step, ": Loglikelihood ", format(estimate$value), 
                ", Parameter: ", private$format_par(estimate$par))
      } else {
        message(" Parameter: ", private$format_par(estimate$par), 
                ", Loglikelihood ", format(old_llh), 
                " -> ", format(estimate$value))
      }
    },
    log_new_rep = function(rep, start_pos) {
      if (private$verbose) {
        message("Repetition ", rep, 
                " starting at ", private$format_par(start_pos), ":")
      }
    },
    log_convergence = function(step) {
      if (private$verbose) message(" Convergence detected")
      private$converged[step] = TRUE
    },
    log_initialization = function(method) {
      if (!private$verbose) return(invisible(NULL))
      message("Determining starting positions using the '", 
              method, "' method")
    },
    log_llh_correction = function() {
      if (private$verbose) message("Correcting likelihoods for best estimates:")
    },
    get_estimates = function(rep) private$estimates[[rep]],
    get_best_estimates = function(n = 5, final = FALSE) {
      "returns the n best estimates per repetition"
      if (final) estimates_list <- list(private$final_estimates)
      else estimates_list <- private$estimates
      best_est <- do.call(rbind, lapply(estimates_list, function(estimates) {
        best_llh <- order(estimates$llh, decreasing = TRUE)[1:n]
        best_llh <- best_llh[!is.na(best_llh)]
        estimates[best_llh, ]
      }))
      best_est[order(best_est$llh, decreasing = TRUE), ]
    },
    create_results = function() {
      "creates the results list the main function returns"
      best_estimate <- self$get_best_estimates(1, TRUE)
      param <- as.numeric(best_estimate[1, -(1:3)])
      res <- list(estimate = private$par_ranges$denormalize(param),
                  loglikelihood = as.numeric(best_estimate[1, 3]),
                  converged = all(private$converged),
                  args = list(repetitions = private$reps,
                              sim = private$sim,
                              max_steps = private$max_steps,
                              init_method = private$init_method))
      class(res) <- c("jaatha_result", class(res))
      res
    }
  )
)

create_jaatha_log <- jaatha_log_class$new
