test.addParameter <- function() {
  dm <- dm.createDemographicModel(11:12, 100)
  dm <- dm.addParameter(dm, "theta", 1, 5)
  dm <- dm.addParameter(dm, "rho", fixed=20)
  checkEquals(dm@parameters$name,  c("theta","rho"))
  checkEquals(dm@parameters$fixed, c(F,T))
  checkEquals(dm@parameters$lower.range, c(1,20))
  checkEquals(dm@parameters$upper.range, c(5,20))
}

test.addFeature <- function() {
  dm <- dm.createDemographicModel(11:12, 100)
  n.feat <- nrow(dm@features)
  dm <- addFeature(dm, type="split", parameter="tau", 
                   lower.range=1, upper.range=10)
  checkEquals(n.feat+1, nrow(dm@features))
  dm <- addFeature(dm, "mutation", "theta", fixed.value=5,
                   pop.source=1, pop.sink=2,
                   time.point="t2", group=3)
  checkEquals(n.feat+2, nrow(dm@features))
}

test.addSummaryStatistics <- function() {
  dm <- dm.createDemographicModel(11:12, 100)
  checkEquals(dm@sum.stats, c('jsfs'))
  dm <- dm.addSummaryStatistic(dm, 'seg.sites')
  checkEquals(dm@sum.stats, c('jsfs', 'seg.sites'))
  checkException(addSummaryStatistic(dm, 'no.existing.sumstat'))
  checkException(addSummaryStatistic(dm, 1:10))
}

test.makeThetaLast <- function() {
  dm <- dm.createDemographicModel(11:12, 100)
  dm <- dm.addMutation(dm,5,20)
  dm <- dm.addSpeciationEvent(dm,0.1,5)
  checkEquals(dm.getParameters(dm)[dm.getNPar(dm)], "theta")
}

test.setMutationModel <- function() {
  dm <- dm.createThetaTauModel(11:12, 100)
  dm <- dm.setMutationModel(dm, "HKY")
  checkTrue("mutation.model" %in% dm@parameters$name)
  checkException(dm <- dm.setMutationModel(dm, "bla"))
}

test.simSumStats <- function() {
  checkException( dm.simSumStats( 1 ) )
  checkException( dm.simSumStats(dm.tt, 1 ) )
  checkException( dm.simSumStats(dm.tt, 1:3 ) )
  checkException( dm.simSumStats(dm.tt, c(2,50)) )

  sum.stats <- dm.simSumStats(dm.tt, c(1,5), "jsfs")
  checkTrue( is.list(sum.stats) )
  checkTrue( !is.null(sum.stats$jsfs) )
  checkTrue( sum(sum.stats$jsfs) > 0 )

  sum.stats <- dm.simSumStats(dm.grp, c(1,5), "jsfs")
  checkTrue( !is.null(sum.stats) )
  checkTrue( !is.null(sum.stats$jsfs.1) )
  checkTrue( sum(sum.stats$jsfs.1) > 0 )
  checkTrue( !is.null(sum.stats$jsfs.2) )
  checkTrue( sum(sum.stats$jsfs.2) > 0 )
  checkTrue( !is.null(sum.stats$jsfs.3) )
  checkTrue( sum(sum.stats$jsfs.3) > 0 )
}

test.parInRange <- function() {
  checkParInRange(dm.tt, c(1, 5))
  checkParInRange(dm.tt, c(2, 7))
  checkParInRange(dm.tt, c(2.1, 7.7))
  checkException( checkParInRange(dm.tt, c(0, 5)) )
  checkException( checkParInRange(dm.tt, c(0, -1)) )
  checkException( checkParInRange(dm.tt, c(10, 1)) )
  checkException( checkParInRange(dm.tt, c(100, 100)) )
  checkException( checkParInRange(dm.tt, 1) )
  checkException( checkParInRange(dm.tt, matrix(1, 2, 2) ))
  checkException( checkParInRange(dm.tt, NULL ))
}

test.addSampleSize <- function() {
  dm <- dm.tt
  checkEquals( 2, sum(dm@features$type == "sample") )
  checkEquals( "11", subset(dm@features, type=="sample" & pop.source==1)$parameter)
  checkEquals( "12", subset(dm@features, type=="sample" & pop.source==2)$parameter)
  checkEquals( 2, nrow(subset(dm@features, type=="sample" & group==0)) )
  checkEquals( 2, nrow(subset(dm@features, type=="sample" & time.point=='0')) )

  dm <- dm.addSampleSize(dm, c(2,5), group=1)
  checkEquals( 4, sum(dm@features$type == "sample") )
  checkEquals( 2, nrow(subset(dm@features, type=="sample" & group==0)) )
  checkEquals( 2, nrow(subset(dm@features, type=="sample" & group==1)) )
}

test.getSampleSize <- function() {
  dm.test <- dm.tt
  checkEquals(c(11L, 12L), dm.getSampleSize(dm.test))

  dm.test <- dm.addSampleSize(dm.test, c(2,5), group=2)
  checkEquals(c(11L, 12L), dm.getSampleSize(dm.test, 0))
  checkEquals(c(11L, 12L), dm.getSampleSize(dm.test, 1))
  checkEquals(c(2L, 5L), dm.getSampleSize(dm.test, 2))
  checkException(dm.getSampleSize(dm.test), silent=TRUE)

  dm.test <- dm.setLociLength(dm.test, 15, group=3)
  checkEquals(c(11L, 12L), dm.getSampleSize(dm.test, 3))
}

test.setLociNumber <- function() {
  dm <- dm.setLociNumber(dm.tt, 17)
  checkEquals( 1, sum(dm@features$type == "loci.number") )
  checkEquals( "17", dm@features$parameter[dm@features$type == "loci.number"] )

  dm <- dm.setLociNumber(dm, 23, group=1)
  checkEquals( 2, sum(dm@features$type == "loci.number") )
  checkEquals( "23",
              dm@features$parameter[dm@features$type=="loci.number"][2] )

  dm <- dm.setLociNumber(dm, 32, group=2)
  checkEquals( 3, sum(dm@features$type == "loci.number") )
}


test.scaleDemographicModel <- function() {
  dm <- dm.setLociNumber(dm.tt, 25, group=1)
  dm <- dm.setLociNumber(dm, 27, group=2)
  dm <- scaleDemographicModel(dm, 5) 
  checkEquals(2L, dm.getLociNumber(dm, 0))
  checkEquals(5L, dm.getLociNumber(dm, 1))
  checkEquals(5L, dm.getLociNumber(dm, 2))
}

test.setLociLength <- function() {
  dm <- dm.setLociLength(dm.tt, 17)
  checkEquals( 1, sum(dm@features$type == "loci.length") )
  checkEquals( "17", dm@features$parameter[dm@features$type == "loci.length"] )

  dm <- dm.setLociLength(dm, 23, group=1)
  checkEquals( 2, sum(dm@features$type == "loci.length") )
  checkEquals( "23",
              dm@features$parameter[dm@features$type=="loci.length"][2] )

  dm <- dm.setLociLength(dm, 32, group=2)
  checkEquals( 3, sum(dm@features$type == "loci.length") )
}

test.getLociNumber <- function() {
  dm.test <- dm.setLociNumber(dm.tt, 17)
  dm.test <- dm.setLociNumber(dm.test, 23, group=1)
  dm.test <- dm.setLociNumber(dm.test, 32, group=2)

  checkEquals(17L, dm.getLociNumber(dm.test, 0))
  checkEquals(23L, dm.getLociNumber(dm.test, 1))
  checkEquals(32L, dm.getLociNumber(dm.test, 2))
  checkException( dm.getLociNumber(dm.test) )

  dm.test <- dm.setLociNumber(dm.tt, 17)
  checkEquals(17L, dm.getLociNumber(dm.test))
  dm.test@features$group[dm.test@features$type == 'loci.number'] <- 1
  checkEquals(17L, dm.getLociNumber(dm.test))
}

test.getLociLength <- function() {
  dm.test <- dm.setLociLength(dm.tt, 17)
  dm.test <- dm.setLociLength(dm.test, 23, group=1)
  dm.test <- dm.setLociLength(dm.test, 32, group=2)

  checkEquals(17L, dm.getLociLength(dm.test, 0))
  checkEquals(23L, dm.getLociLength(dm.test, 1))
  checkEquals(32L, dm.getLociLength(dm.test, 2))
  checkException( dm.getLociLength(dm.test) )

  dm.test <- dm.setLociLength(dm.tt, 17)
  checkEquals(17L, dm.getLociLength(dm.test))
  dm.test@features$group[dm.test@features$type == 'loci.length'] <- 1
  checkEquals(17L, dm.getLociLength(dm.test))
}

test.getGroups <- function() {
  checkEquals(1, dm.getGroups(dm.tt))
  dm <- dm.setLociLength(dm.tt, 23, group=1)
  checkEquals(1, dm.getGroups(dm))
  dm <- dm.setLociLength(dm.tt, 32, group=2)
  checkEquals(1:2, dm.getGroups(dm))
  dm <- dm.setLociNumber(dm, 32, group=4)
  checkEquals(c(1:2,4), dm.getGroups(dm))
}

test.searchFeature <- function() {
  checkEquals(2, nrow(searchFeature(dm.tt, type='sample')))
  checkEquals(3, nrow(searchFeature(dm.tt, type=c('split', 'sample'))))
  checkEquals(2, nrow(searchFeature(dm.tt, type=c('split', 'sample'),
                                    pop.sink=NA)))
  checkEquals(1, nrow(searchFeature(dm.tt, time.point='tau', group=0)))
  checkEquals(3, nrow(searchFeature(dm.tt, time.point=NA)))
}

test.generateGroupModel <- function() {
  dm <- generateGroupModel(dm.tt, 1)
  checkEquals(dm.tt@features, dm@features)

  dm <- dm.setLociLength(dm.tt, 23, group=1)
  dm <- generateGroupModel(dm, 1)
  checkEquals(nrow(dm.tt@features), nrow(dm@features))
  checkEquals(1, sum(dm@features$group == 1))
  checkEquals('23', dm@features$parameter[dm@features$group == 1])

  dm.3 <- dm.setLociLength(dm.tt, 23, group=1)
  dm.3 <- dm.setLociLength(dm.3, 30, group=2)
  dm.3 <- dm.setLociNumber(dm.3, 31, group=2)
  dm <- generateGroupModel(dm.3, 1)
  checkEquals(nrow(dm.tt@features), nrow(dm@features))
  checkEquals(1, sum(dm@features$group == 1))
  checkEquals('23', dm@features$parameter[dm@features$group == 1])
  dm <- generateGroupModel(dm.3, 2)
  checkEquals(nrow(dm.tt@features), nrow(dm@features))
  checkEquals(0, sum(dm@features$group == 1))
  checkEquals(2, sum(dm@features$group == 2))
}

test.finalizeDM <- function() {
  finalizeDM(dm.tt)
  finalizeDM(dm.hky)
  finalizeDM(dm.grp)
}

## -- Fixed bugs ----------------------------------------
test.showEmptyModel <- function() {
  dm <- dm.createDemographicModel(25:26, 100)
  invisible(print(dm))
}
