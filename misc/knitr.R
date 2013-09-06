#!/usr/bin/Rscript --vanilla
#
# knitr.R
# Script to call knitr on a *.Rnw file with a specific cache dir
#
# Usage: knitr.R <file-to-knit> [<cache-dir>]
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2013-09-06
# Licence:  GPLv3 or later
#

cmd.args <- commandArgs(TRUE)

library('knitr')

if (!is.na(cmd.args[2])) {
  if (!file.exists(cmd.args[2])) stop("cache dir not found: ", cmd.args[2])
  opts_chunk$set(cache.path=cmd.args[2])
} 

if (is.na(cmd.args[1]) | !file.exists(cmd.args[1])) stop("File to knit not found: ",
                                                 cmd.args[1])
knit2pdf(cmd.args[1], encoding='UTF-8')
