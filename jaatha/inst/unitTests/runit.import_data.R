#!/usr/bin/Rscript --vanilla
#
# ../inst/unitTests/runit.import_data
# %DESCRIPTION%
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-08-28
# Licence:  GPLv3 or later
#

library(RUnit)

marker.table <- rbind(c("A","A","G","A","A","A","A"),  #Valid Snp
                      c("A","A","G","A","A","G","G"),  #Valid Snp
                      c("A","A","G","G","A","G","T"),  #Valid Snp
                      c("A","G","A","A","A","G","G"),  #Outg disagree
                      c("A","A","A","A","A","A","A"))  #fixed

outg.mask <- c(T, T, F, F, F, F, F)
pop1.mask <- c(F, F, T, T, T, F, F)
pop2.mask <- c(F, F, F, F, F, T, T)

test.isSNP <- function(){
  snp.mask <- apply(marker.table, 1, isSnp, outg.mask=outg.mask)
  checkEquals(snp.mask, c(T, T, T, F, F))
}

test.getSnpTypes <- function() {
  snp.types <- getSNPTypes(marker.table[1:3, ], pop1.mask, pop2.mask, folded=F)
  checkEquals(snp.types[1, ], c(1, 0))
  checkEquals(snp.types[2, ], c(1, 2))
  checkEquals(snp.types[3, ], c(2, 2))
}

test.readMarkerTable <- function(){
  jsfs <- markerTableToJsfs(marker.table, 3:5, 6:7, 1:2)
  checkEquals(jsfs[2,1], 1)
  checkEquals(jsfs[2,3], 1)
  checkEquals(jsfs[3,3], 1)
  checkEquals(sum(jsfs), 3)
}
