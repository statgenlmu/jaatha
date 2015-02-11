# Example for profiling compiled code
# Needs googe-perftools and littler installed
# Calls with
# LD_PRELOAD="/usr/lib/libprofiler.so.0" CPUPROFILE=/tmp/rprof.log r test.R
# and the view with e.g.
# google-pprof --gv /usr/bin/r /tmp/rprof.log
# or
# google-pprof --cum --text /usr/bin/r /tmp/rprof.log | less

# First create needed data and the comment out again
#
# dm <- dm.createThetaTauModel(c(20,25), 50000)
# sim_data <- dm.simSumStats(dm.addSummaryStatistic(dm, 'seg.sites'), c(.1,10))
# llm <- dm.getLociLengthMatrix(dm)
# save(llm, sim_data, file='misc//test_data.Rda')

load('test_data.Rda')
jaatha:::calcPercentFpcViolation(sim_data$seg.sites, 1:20, llm)
