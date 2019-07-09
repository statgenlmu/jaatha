jaatha 3.2.1
============

* Fixes a problem with a vector-valued condition during error handling which
  became appartent in the tests.
* Exports `coarsen_jsfs` function (#141).



jaatha 3.2.0
============

* Adds the `cores_per_run` to `boot_jaatha` which can be set to use more than 
  one CPU core for each replica (#128).
* Coala support: Adds JSFS summary options `folded_sums`, which calculates
  the standard JSFS coarsening from an unpolarized spectrum (#130).
* Now tries to continue the estimation even when a few simulations return 
  errors (#131).
* Skips tests if `testthat` is not available (#132).
* Adds new options `final_sim` and `zoom_in_steps` to the main jaatha `method` 
  (#134).
* Fixes a bug that lead to a broken assertion when jaatha was used on a model
  with one parameter and `repetitions` was greater than one (#135). 
  Thanks to Amaya Romero for reporting and locating this bug!.



jaatha 3.1.1
============

* Fixes broken suppression of messages with `verbose` argument of the `jaatha` 
  function (#118).
* Removes a buggy check for summary statistics that are insensitve about the 
  parameters introduced in version 3.1.0. Together with the previous bug, this
  could lead to false results in some applications of _jaatha_. For any runs
  conducted with that version, check if the output contained an estimated 
  likelihood of `1` and repreat the runs with the new version if so (#124).
* Adds more assertions to the code (#120).



jaatha 3.1.0
============

* Adds a limit to the number of cached simulations (#103).
* Adds the `random` initialization method for using random starting positions (#104).
* Also captures `model` and `data` in the `pars` argument of an jaatha analysis
  and removed the corresponding arguments from `boot_jaatha` (#105).
* Adds a custom print methods for Jaatha's results (#106).
* Optimizes scheduling of simulations based on their duration. If simulations
  take more than five seconds, the `mc.preschedule` argument of 
  `parallel::mclapply` is now set to `FALSE` (#107).
* Adds the `block_width` argument to `jaatha` (#108).
* Also exports `create_jaatha_model.function` and  
  `create_jaatha_model.coalmodel` as functions (#109).
* Increases robustsness to not-converging GLMs. If a GLM fails to converge, now
  only the current repetition is aborded instead of the compelete analysis 
  (#110, #111). 
* Dump frames to tempdir() on simulation errors (#113).



jaatha 3.0.0
============

* This version is a complete rewrite of the package. The R interface completely
  redesigned, and should now be more similar to other R packages.
* All functionality for coalescent simulation was removed from the package,
  and is now available via the separate package `coala`.
* The use of `jaatha` as a general method is now more prominent.
* The documentation was revised. The vignettes are now included in the package.
* More methods for generating starting positions were added to the package.
* Jaatha now uses the `boot` package for bootstrapping calculations.
* Internally, all S4 classes where replaced with R6 classes to increase the
  performance.



jaatha 2.7.0
============

* Implement a candidate loci approach for selection. So far only for internal
  testing.
* Allow for gamma distributed recombination and mutation rates
  between loci for genomic parameters (mutation and recombination rates).
* Bug fix: Fix two bugs with calculation of confidence intervals.
* Fix a harmless memory leak in experimental code (#19)



jaatha 2.6.0
============

* Switched unittests from RUnit to testthat
* Added a Rcpp function to parse tree from ms(ms) output for seq-gen. This should
  make seq-gen less platform dependent (it still requires a Unix platform).
* Add a subset parameter to Jaatha.confidenceIntervals that allows to distribute
  the calculation on multiple machines.
* Bug fix: in 2.5.1, the ms output file was not delete when using seq-gen. This can cause
  Jaatha to fail if the disk runs out of space (#14). 
* Bug fix: Fixed a severe bug in the `Jaatha.confidenceIntevals` function. There often
  parameters different from the maximum likelihood estimates were used to
  calculate the confidence intervals (#15). 



jaatha 2.5.1
============

* Made it possible to use seq-gen together with msms (Linux only).
* Bug fix: Fixed a bug that caused seq-gen to fail if the mutation parameter was not the
  last model parameter.



jaatha 2.5.0
============

* Fixed wrong description of N0 as size of the ancestral population in the documentation (#4)
* Fixed an bug with made jaatha abort when a point was slightly outside
  parameter range (#5)
* Fixed bug in calculation of JSFS introduced in 2.4-1 (#7)
* Now old simulations are actually reused (#8)
* It is no longer necessary to have a simulation program that can simulate all
  groups (#6)
* Different groups now really use different simulation programs
* Simplified internal simulation program interface
* Print dispatch for demographic model
* Modularized parsing of simulation output
* Models in different groups can now use different simulation software
* Now Jaatha simulates more data if it fails to fit the GLM
* Allow GLM more steps to converge
* Different Groups can now have different summary statistics
* Prioritizing system for simulations programs
* Removed additional tmp-folders
* Removed dependency on foreach
* Changed dependency from dplyr to reshape2
* Removed optional logging into a file



jaatha 2.4.0
============

* Jaatha.initialSearch() now divides the parameter space by all axis in turn,
  which makes the number of block linear instead of exponential in the
  parameter number.  
* Add an optional border to smoothing type summary statistics with is fitted
  via traditional transformations.
* Using `dplyr` instead of `plyr` for converting JSFS into data frames because
  it is significantly faster.
* Striped unused components from the `glm` objects.



jaatha 2.3.0
============

* Removed `epsilon`. We now stop the refined search after the
  likelihood did not improve for 10 consecutive steps.
* Removed `weight` parameter to down-weight old simulations,
  because it had no measurable effect on the estimates.
* Support for multiple groups of loci
* Possibility to smooth arrays of arbitrary dimension



jaatha 2.2.0
============

* Added function for calculating BCa confidence intervals
* Added `rerun` option to initial & refined search 
* Switched form doMC to multicore for parallelization
* Support for multiple independent summary statistics, even though this is not 
  yet used in the default front end.
* Differentiated between unit and integration tests and moved the unit tests
  back into the package.
* Reworked the reusing of simulations from previous steps in the refined search.
  Now only simulations that actually are inside the current block are reused. 
* Removed the "score". We are always using the log-likelihood now.
* Added support for a new "poisson.smoothing" type of summary statistics. 
* Fixed problem with naming of temporary directories if Jaatha was executed
  multiple times in different threads
* Fixed problem with Rcpp in size_t on 64 bit windows
* Fixed warning messages in cpp files
* Updated Makevars



jaatha 2.1.0
============

* Removed the "externalTheta" option to simplify further development
* Rewrote and cleaned up many of the internal functions
* Further separated Jaatha and Demographic Model to simplify use of Jaatha
  in other areas
* Added caching for log factorials
* Ported msFile2jsfs to Rcpp
* Removed vignettes from the package due to their long build time; The pdfs will
  be provided on Jaatha's homepage instead 
* Added citation information to the package  
* Added link to github repo to package description
* Removed unused optimization strategies
* Bugfix: Changed the example fasta file for compatibility with next ape version
* Bugfix: Fixed encoding issue with bibtex for the packages vignette
* Bugfix: Removed debug functions in C++ files that prevented installation on R-2.14



jaatha 2.0.0
============

* Updated Documentation and Vignette
* New Feature: Added caching of ms and seqgen commands
* Changed Feature: Made use of /dev/shm optional
* Bugfix: Multiple tmp-dirs used during multicore refined search
* Bugfix: Multiple searches for seqgen binary if running on multiple cores
* Fixed a bug with one parameter models
* Added support for seqgen
* GLM fitting algorithm can now take more steps to converge
* Updated the vignette
* Ported seqgenFile2Jsfs.cc to Rcpp
* Added parallelization using doMC.
* Some improvements to the UI.
* Updated copyright information.
* Added import from ape.
* Updated the vignette.
* Use /dev/shm for temp files if available
* Added "folded" option to use a folded JSFS
* Added documentation for fasta2jsfs()
* Made fasta2jsfs calculate the folded JSFS if no outgroup is given.
* Reworked the demographic model class to be more flexible and consistent.

