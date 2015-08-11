jaatha_log_class <- R6Class("jaatha_log", 
  private = list(
    estimates = NULL,
    final_estimates = NULL,
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
      private$estimates <- lapply(1:reps, function(i) {
        estimates <- matrix(NA, max_steps, par_number + 3)
        colnames(estimates) <- c("rep", "step", "llh", par_names)
        as.data.frame(estimates)
      })
      private$final_estimates <- private$estimates[[1]][rep(1, 5 * par_number), ]
      private$max_steps <- max_steps
      private$verbose <- verbose
      private$par_ranges <- model$get_par_ranges()
    },
    log_estimate = function(rep, step, estimate) {
      private$estimates[[rep]][step, ] <- 
        c(rep, step, estimate$value, estimate$par)
      
      if (private$verbose) {
        message("Step ", step, ": Loglikelihood ", format(estimate$value), 
                ", Parameter: ", private$format_par(estimate$par))
      }
    },
    log_final_estimate = function(i, prev, llh) {
      prev$llh <- llh
      private$final_estimates[i, ] <- prev
        
      if (private$verbose) {
        message("Final llh correction", i, ": Loglikelihood ", format(llh), 
                ", Parameter: ", private$format_par(prev[1, -(1:3)]))
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
    create_results = function() {}
  )
)

create_jaatha_log <- jaatha_log_class$new
