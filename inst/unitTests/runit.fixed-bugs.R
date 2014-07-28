# unit_tests/runit.fixed-bugs
# Unit test ensuring regressions of already fixed bugs
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-11-30
# Licence:  GPLv3 or later
#

## -- Fixed bugs ----------------------------------------

#print() failed for empty demographic models
test.printEmptyDM <- function() {
  dm <- dm.createDemographicModel(25:26, 100)
  print(dm)
}

test.printGroupDM <- function() {
  print(dm.grp)
}

test.ApeCodeInVignette <- function() {
  require('ape') || stop('Package "ape" is required for unit tests')
  sample.file <- system.file('example_fasta_files/sample.fasta', package='jaatha')
  sample.data <- read.dna(sample.file, format='fasta', as.character=TRUE)
  sample.jsfs <- calculateJsfs(sample.data,
                               pop1.rows=3:7, 
                               pop2.rows=8:12,
                               outgroup.row=1:2)
  checkTrue(sum(sample.jsfs) > 0)
}
