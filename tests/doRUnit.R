#!/usr/bin/Rscript --vanilla
#
# Script to run Jaatha's unit tests
# 
# Taken form:
# http://rwiki.sciviews.org/doku.php?id=developers:runit
#

stopifnot(require("RUnit", quietly=TRUE))

run.integration.tests <- Sys.getenv("INTEGRATION_TESTS") == "TRUE"
rcmdcheck <- Sys.getenv("RCMDCHECK") != "FALSE"
pkg <- "jaatha"



## --- Setup ---

if(rcmdcheck) {
  path.unit.tests <- system.file(package=pkg, "unitTests")
  path.integration.tests <- system.file(package=pkg, "integrationTests")
} else {
  path.unit.tests <- file.path(getwd(), "inst", "unitTests")
  path.integration.tests <- file.path(getwd(), "inst", "integrationTests")
}


library(package=pkg, character.only=TRUE)

## load the name space to allow testing of private functions
if (is.element(pkg, loadedNamespaces()))
  attach(loadNamespace(pkg), name=paste("namespace", pkg, sep=":"), pos=3)

test.setup <- paste(path.unit.tests, "test-setup.Rda", sep="/")
if (!file.exists(test.setup)) stop("Failed to load ", test.setup) 
load(test.setup)
rm(test.setup)

## ---  Set Options ---
runit.options <- getOption("RUnit")
runit.options$silent <- TRUE
options(list(RUnit=runit.options))

## ---  Define tests ---
test.suites <- list()
test.suites[['unit.tests']] <- defineTestSuite(name="Unit Tests",
                                               dirs=path.unit.tests)

if (run.integration.tests) {
  test.suites[['integration.tests']] <- defineTestSuite(name="Integration Tests",
                                                        testFileRegexp = "^it.+\\.[rR]$",
                                                        dirs=path.integration.tests)
}



## --- Run tests ---
tests <- runTestSuite(test.suites)



## --- Evaluate Results ---
cat("------------------- TEST SUMMARY ---------------------\n\n")
printTextProtocol(tests, showDetails=FALSE)
tmp <- getErrors(tests)

if(tmp$nFail > 0 | tmp$nErr > 0) {
  stop(paste("\n\nunit testing (#test failures: ", tmp$nFail,
             ", #R errors: ",  tmp$nErr, ")\n\n", sep=""))
}
