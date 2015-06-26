Stat_sfs <- R6Class("Stat_SFS", 
  inherit = Stat_PoiInd,
    public = list(
      initialize = function(seg_sites, model, stat) {
        name <- stat$get_name()
        sfs <- stat$calculate(seg_sites, NULL, NULL, model)
        fake_sim_data <- list()
        fake_sim_data[[name]] <- sfs
        super$initialize(fake_sim_data, stat$get_name())
      },
      transform = function(sim_data) sim_data[[private$name]]
   )
)