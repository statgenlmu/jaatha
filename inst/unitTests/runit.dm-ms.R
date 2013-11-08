dm <- dm.createThetaTauModel(c(10,25), 100)

test.callMS <- function() {
  checkException(callMs())
  checkTrue(file.exists(callMs("-t 5", dm)))
}

test.msSingleSimFunc <- function() {
  sum.stats <- msSingleSimFunc(dm, c(1,10))
  checkTrue(is.list(sum.stats))
  checkTrue(is.array(sum.stats$jsfs))
  checkTrue(sum(sum.stats$jsfs) > 0)
}

test.simProg <- function(){
  checkTrue(!is.null(.jaatha$simProgs[["ms"]]))
}
