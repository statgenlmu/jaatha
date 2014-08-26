#!/usr/bin/Rscript
library(testthat)
library(jaatha)

if (file.exists("DESCRIPTION")) setwd("tests//integration")
source("..//testthat//test-aaa-setup.R")
test_dir('.', reporter='stop')