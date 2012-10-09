#!/usr/bin/Rscript --vanilla
#
# doRUnit.R
# Script to run Jaathas Unit tests
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-10-08
# Licence:  GPLv3 or later
#

# Load required packages
require("RUnit", quietly=TRUE)
require("devtools", quietly=TRUE)

# Install the current version of Jaatha
#install("../jaatha", quick=T)
require("jaatha")

# Load Package Namespace
if (is.element("jaatha", loadedNamespaces()))
  attach(loadNamespace("jaatha"), name=paste("namespace", "jaatha", sep=":"), pos=3)

## --- Testing ---

## Define tests
testSuite <- defineTestSuite(name="Jaatha Unit Testing", dirs=".")

## Run
tests <- runTestSuite(testSuite)

## Default report name
pathReport <- file.path("./results/")
if (!file.exists(pathReport)) dir.create(pathReport)

## Report to stdout and text files
cat("------------------- UNIT TEST SUMMARY ---------------------\n\n")
printTextProtocol(tests, showDetails=FALSE)
printTextProtocol(tests, showDetails=FALSE,
                  fileName=paste(pathReport, "summary.txt", sep=""))
printTextProtocol(tests, showDetails=TRUE,
                  fileName=paste(pathReport, "detailed.txt", sep=""))

## Report to HTML file
printHTMLProtocol(tests, fileName=paste(pathReport, "detailed.html", sep=""))

## Return stop() to cause R CMD check stop in case of
##  - failures i.e. FALSE to unit tests or
##  - errors i.e. R errors
tmp <- getErrors(tests)
if(tmp$nFail > 0 | tmp$nErr > 0) {
  stop(paste("\n\nunit testing failed (#test failures: ", tmp$nFail,
             ", #R errors: ",  tmp$nErr, ")\n\n", sep=""))
}
