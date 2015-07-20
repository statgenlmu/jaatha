do_initial_search <- function(model, data, repetititons, sim, cores) {
  
  par_number <- model$get_par_ranges()$get_par_number()
  blocks_per_par <- determine_bpp(par_number, repetititons)
  blocks <- create_initial_blocks(model$get_par_ranges(), blocks_per_par)
  
  
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
    return(list(block_class$new(basic_block)))
  }
  
  blocks <- list()
  length(blocks) <- par_number * blocks_per_par
  for (i in 1:par_number) {
    for (j in 1:blocks_per_par) {
      new_block <- basic_block
      new_block[i, 2] <- 1 / blocks_per_par
      new_block[i, ] <- new_block[i, ] + (j - 1) / blocks_per_par
      blocks[[(j - 1) * par_number + i]] <- block_class$new(new_block)
    }
  }
  
  blocks
}
