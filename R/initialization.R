get_start_pos <- function(model, data, reps, sim, init_method, cores, 
                          sim_cache) {
  
  start_pos <- NULL
  if (init_method[1] == "zoom-in") {
    start_pos <- do_zoom_in_search(model, data, reps, sim, cores, sim_cache)
  } else if (init_method[1] == "initial-search") {
    start_pos <- do_initial_search(model, data, reps, sim, cores, sim_cache)
  } else if (init_method[1] == "middle") {
    start_pos <- matrix(.5, reps, model$get_par_number())
  } else {
    stop("Unknown initialization method: ", init_method[1])
  }
  
  assert_that(is.matrix(start_pos))
  assert_that(all((dim(start_pos) == c(reps, model$get_par_number()))))
  start_pos
}


do_initial_search <- function(model, data, reps, sim, cores, sim_cache) {
  # Divide the parameter space in blocks
  par_number <- model$get_par_ranges()$get_par_number()
  blocks_per_par <- determine_bpp(par_number, reps)
  blocks <- create_initial_blocks(model$get_par_ranges(), blocks_per_par)
  
  # Get an estimate infor each block
  estimates <- lapply(blocks, estimate_local_ml, model, data, 
                      sim, cores, sim_cache)
  
  # Return the parameters for the best estimates
  best_indexes <- order(vapply(estimates, function(x) x$value, numeric(1)), 
                        decreasing = TRUE)[1:reps]
  t(vapply(estimates[best_indexes], function(x) x$par, 
           numeric(model$get_par_number())))
}


determine_bpp <- function(par_number, repetitions) {
  "Finds the minimal number of blocks per par (bpp) needed"
  if (repetitions == 1) return(1)
  blocks_per_par <- 2
  while (repetitions > blocks_per_par * par_number) {
    blocks_per_par <- blocks_per_par + 1
  }
  blocks_per_par
}


create_initial_blocks <- function(par_ranges, blocks_per_par) {
  par_number <- par_ranges$get_par_number()
  basic_block <- matrix(c(0, 1), par_number, 2, byrow = TRUE)

  if (blocks_per_par == 1) {
    return(list(create_block(basic_block)))
  }
  
  blocks <- list()
  length(blocks) <- par_number * blocks_per_par
  for (i in 1:par_number) {
    for (j in 1:blocks_per_par) {
      new_block <- basic_block
      new_block[i, 2] <- 1 / blocks_per_par
      new_block[i, ] <- new_block[i, ] + (j - 1) / blocks_per_par
      blocks[[(j - 1) * par_number + i]] <- create_block(new_block)
    }
  }
  
  blocks
}


do_zoom_in_search <- function(model, data, reps, sim, cores, sim_cache) {
  t(vapply(1:reps, function(i) {
    middle <- rep(.5, model$get_par_number())
    for (block_width in c(1, 0.5, 0.25)) {
      block <- create_block(cbind(middle - block_width * .5,
                                  middle + block_width * .5), cut = TRUE)
      middle <- estimate_local_ml(block, model, data, sim, cores, sim_cache)$par
    }
    middle
  }, numeric(model$get_par_number())))
}
