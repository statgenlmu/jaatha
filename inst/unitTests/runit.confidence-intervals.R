test.convertSimDataToSumStats <- function() {
  sim.data <- sim.data.csi[[1]]
  sim.data$poisson.vector <- rep(7.5, 6)
  sum.stats <- convertSimDataToSumStats(sim.data, jaatha.csi@sum.stats)

  checkTrue( length(sum.stats) == length( jaatha.csi@sum.stats ) )
  checkTrue( sum.stats$poisson.vector$method == 
             jaatha.csi@sum.stats$poisson.vector$method )
  checkTrue( all(sum.stats$poisson.vector$value == 7.5) )
  checkTrue( all(sum.stats$poisson.vector$value.transformed == 7.5) )
}
