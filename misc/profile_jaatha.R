# Example script for profiling the initiual and refined search.

model <- create_test_model()
data <- create_test_data(model)

proffile <- tempfile('jaatha_is2.prof')

gc()
set.seed(17)

Rprof(proffile, memory.profiling = TRUE)
results <- jaatha(model, data, repetitions = 2, sim = 500, cores = 2)
Rprof(NULL)

summaryRprof(filename = proffile, memory = 'both')



# devtools::install_github("hadley/lineprof")
library(lineprof)
pr <- lineprof(results <- jaatha(model, data, repetitions = 2, sim = 100, cores = 1))
shine(pr)
