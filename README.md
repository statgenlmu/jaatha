[![Linux Build Status](https://travis-ci.org/statgenlmu/jaatha.svg?branch=master)](https://travis-ci.org/statgenlmu/jaatha) 
[![Windows Build status](https://ci.appveyor.com/api/projects/status/g4adpum1pkyn7ajn/branch/master?svg=true)](https://ci.appveyor.com/project/paulstaab/jaatha/branch/master)
[![Coverage Status](https://coveralls.io/repos/statgenlmu/jaatha/badge.svg?branch=master&service=github)](https://coveralls.io/github/statgenlmu/jaatha?branch=master)
[![CRAN Status](http://www.r-pkg.org/badges/version/jaatha)](https://cran.r-project.org/web/packages/jaatha/)


Jaatha
======

Jaatha is an estimation method that uses computer simulations to produce
maximum-likelihood estimates even when the likelihood function can not be
evaluated directly. It can be applied whenever it is feasible to conduct many
simulations, but works best when the data is at least approximately Poisson
distributed.

Jaatha was originally designed for demographic inference in evolutionary
biology. It has optional support for conducting coalescent simulation using
the [coala](https://github.com/statgenlmu/coala) R package, but can also be 
used for different applications.

Jaatha is implemented as an [R](https://www.r-project.org) package and available on
[CRAN](https://cran.r-project.org/web/packages/jaatha/).



Installation
------------

Jaatha can be installed from CRAN using the `install.packages` command:

```R
install.packages('jaatha')
```


Usage
-----

The R package includes an introduction vignette that explains how a jaatha
analysis is conducted. After the package is installed, you can open the
vignette using:

```R
vignette("jaatha-intro")
```

A second vignette called `jaatha-evolution` describes how jaatha can be used
together with `coala` for demographic inference.

Further help is provided using R's help system, in particular via `?jaatha`,
`?create_jaatha_model` and `?create_jaatha_data`.



References
----------

Jaatha's original algorithm is described in the publication:

> L. Naduvilezhath, L.E. Rose and D. Metzler:
> Jaatha: a fast composite-likelihood approach to estimate demographic 
> parameters. Molecular Ecology 20(13):2709-23 (2011).

The revised version of the algorithm that is implemented in this package 
is described in:

> L.A. Mathew, P.R. Staab, L.E. Rose and D. Metzler:
> [Why to account for finite sites in population genetic studies and 
> how to do this with Jaatha 2.0](http://onlinelibrary.wiley.com/doi/10.1002/ece3.722/abstract). 
> Ecology and Evolution (2013).



Development
-----------

Jaatha is developed openly on [GitHub](https://github.com/statgenlmu/jaatha). 
Feel free to open an issue there if you encounter problems using Jaatha or 
have suggestions for future versions.

The current development version can be installed using:

```R
devtools::install_github('statgenlmu/jaatha')
```
