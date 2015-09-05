[![Linux Build Status](https://travis-ci.org/statgenlmu/jaatha.svg?branch=master)](https://travis-ci.org/statgenlmu/jaatha) 
[![Windows Build status](https://ci.appveyor.com/api/projects/status/g4adpum1pkyn7ajn/branch/master?svg=true)](https://ci.appveyor.com/project/paulstaab/jaatha/branch/master)
[![Coverage Status](https://coveralls.io/repos/statgenlmu/jaatha/badge.svg?branch=master&service=github)](https://coveralls.io/github/statgenlmu/jaatha?branch=master)
[![CRAN Status](http://www.r-pkg.org/badges/version/jaatha)](http://cran.r-project.org/web/packages/jaatha)


Jaatha
======

Jaatha is a composite-likelihood, simulation-based parameter estimation method. 
It is designed for situations where the data are assumed to be a vector of independent 
and approximately Poisson distributed variables, but the expectations of the distributions 
depend on the model parameters in an unknown way. Hence, it is not possible to 
evaluate the likelihood function directly. Jaatha fits local Generalized Linear
Models to simulated data to approximate the (composite-) likelihood function and
to produce maximum likelihood estimates.

The method was originally developed to infer the evolutionary history of
two populations from genetic data, but is not limited to this settings. Practical 
instructions on running Jaatha are provided in the `Jaatha-Introduction` vignette. 

Jaatha is developed openly on [GitHub][https://github.com/statgenlmu/jaatha]. 
Feel free to open an issue there if you encounter problems using Jaatha or have 
suggestions for future versions.


Installation
------------

### Stable Version

To install the current stable version of jaatha from CRAN, type

```R
install.packages('jaatha')
```

in R.


### Development Version

You can install the development version from [GitHub][4] using: 

```R
devtools::install_github('statgenlmu/jaatha')
```



Citation
--------
If you use Jaatha in a scientific publication, please cite

> L.A. Mathew, P.R. Staab, L.E. Rose and D. Metzler:
> [Why to account for finite sites in population genetic studies and 
> how to do this with Jaatha 2.0][http://onlinelibrary.wiley.com/doi/10.1002/ece3.722/abstract]. 
> Ecology and Evolution (2013).



Links
-----

* [Jaatha's Homepage](http://evol.bio.lmu.de/_statgen/software/jaatha)
* [Source Code on GitHub](https://github.com/statgenlmu/jaatha)
* [Bug tracker](https://github.com/paulstaab/statgenlmu/issues)
* [Jaatha's page on CRAN](http://cran.r-project.org/web/packages/jaatha/index.html)

