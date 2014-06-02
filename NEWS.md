Version 2.4
==================
## Algorithm changes:
- `Jaatha.initialSearch()` now divides the parameter space by all axis in turn,
  which makes the number of block linear instead of exponential in the
  parameter number.  

## New Features:
- Add an optional border to smoothing type summary statistics with is fitted
  via traditional transformations.

## Improvements:
- `dplyr` instead of `dplyr` for converting JSFS into data frames as the latter 
  is significantly faster.
- Striped unused components from the `glm` objects.


Version 2.3
==================
## Algorithm changes: 
  - Removed `epsilon`. We now stop the refined search after the
    likelihood did not improve for 10 consecutive steps.
  - Removed `weight` parameter to down-weight old simulations,
    because it had no measurable effect on the estimates.

## New Features: 
  - Support for multiple groups of loci
  - Possibility to smooth arrays of arbitrary dimension

## Improvements: 
  - Clarifies 'The Jaatha HowTo'


Version 2.2
==================
+ Added function for calculating BCa confidence intervals
+ Added `rerun` option to initial & refined search 
+ Switched form doMC to multicore for parallelization
+ Support for multiple independent summary statistics, even though this is not 
  yet used in the default front end.
+ Differentiated between unit and integration tests and moved the unit tests
  back into the package.
+ Reworked the reusing of simulations from previous steps in the refined search.
  Now only simulations that actually are inside the current block are reused. 
+ Removed the "score". We are always using the log-likelihood now.
+ Added support for a new "poisson.smoothing" type of summary statistics. 
+ Fixed problem with naming of temporary directories if Jaatha was executed
multiple times in different threads
+ Fixed problem with Rcpp in size_t on 64 bit windows
+ Fixed warning messages in cpp files
+ Updated Makevars


Version 2.1
==================
+ Removed the "externalTheta" option to simplify further development
+ Reworte and cleaned up many of the internal functions
+ Further seperated Jaatha and Demographic Model to simplify use of Jaatha
in other areas
+ Added caching for log factorials
+ Ported msFile2jsfs to Rcpp
+ Removed vignettes from the package due to their long build time; The pdfs will
be provided on Jaatha's homepage instead 
+ Added citation information to the package  
+ Added link to github repo to package description
+ Removed unused optimization strategies
+ Bugfix: Changed the example fasta file for compatibility with next ape version
+ Bugfix: Fixed encoding issue with bibtex for the packages vignette
+ Bugfix: Removed debug functions in C++ files that prevented installation on R-2.14


Version 2.0.0
==================
+ Updated Documentation and Vignette
+ New Feature: Added caching of ms and seqgen commands
+ Changed Feature: Made use of /dev/shm optional
+ Bugfix: Multiple tmp-dirs used during multicore refined search
+ Bugfix: Multiple searches for seqgen binary if running on multiple cores
+ Fixed a bug with one parameter models
+ Added support for seqgen
+ GLM fitting algorithm can now take more steps to converge
+ Updated the vignette
+ Ported seqgenFile2Jsfs.cc to Rcpp
+ Added parallelization using doMC.
+ Some improvements to the UI.
+ Updated copyright information.
+ Added import from ape.
+ Updated the vignette.
+ Use /dev/shm for temp files if available
+ Added "folded" option to use a folded JSFS
+ Added documentation for fasta2jsfs()
+ Made fasta2jsfs calculate the folded JSFS if no outgroup is given.
+ Reworked the demographic model class to be more flexible and consistent.
