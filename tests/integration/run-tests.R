#!/usr/bin/Rscript
library(testthat)
library(devtools)
library(methods)

# Install jaatha & attach its namespace
if (file.exists("DESCRIPTION")) setwd("tests//integration")
install.packages('../..', repos=NULL)
attach(loadNamespace('jaatha'), name=paste("namespace", 'jaatha', sep=":"), pos=3)

source("..//testthat//test-aaa-setup.R")
test_dir('.', reporter='stop')
