Jaatha's README
===============

Jaatha is a fast parameter estimation method primarily designed for Evolutionary
Biology. The method is described in 

> Mathew, L. A., Staab, P. R., Rose, L. E. and Metzler, D. (2013), 
> [Why to account for finite sites in population genetic studies and 
> how to do this with Jaatha 2.0][1]. Ecology and Evolution.

Practical instructions for running Jaatha are provided in the 
[The Jaatha HowTo][2]. Instructions how to use Jaatha with a non-standard 
simulation method are given in the [Custom Simulation Method HowTo][3]. 

Jaatha is developed openly on [GitHub][4]. Feel free to open an issue there if 
you encounter problems using Jaatha or have suggestions for future versions.


Installation
------------

## Stable Version from CRAN:

[![Build Status](https://travis-ci.org/paulstaab/jaatha.png?branch=master)](https://travis-ci.org/paulstaab/jaatha)

To install the current stable version of jaatha from CRAN, simply type

```R
install.packages('jaatha')
```

in R.

## Development Version from GitHub  

[![Build Status](https://travis-ci.org/paulstaab/jaatha.png?branch=develop)](https://travis-ci.org/paulstaab/jaatha)

You can install the development version from GitHub using: 

```R
devtools::install_github('paulstaab/jaatha', ref='develop')
```


Usage
-----

Please refer to the [The Jaatha HowTo][2] for usage information.


Links
-----

[1]: http://onlinelibrary.wiley.com/doi/10.1002/ece3.722/abstract
[2]: https://github.com/paulstaab/jaatha/raw/master/howtos/jaatha_howto.pdf
[3]: https://github.com/paulstaab/jaatha/raw/master/howtos/custom_simulator_howto.pdf
[4]: https://github.com/paulstaab/jaatha

* [Jaatha's Homepage](http://evol.bio.lmu.de/_statgen/software/jaatha)
* [Source Code on GitHub](https://github.com/paulstaab/jaatha)
* [Bug tracker](https://github.com/paulstaab/jaatha/issues)
* [Jaatha's page on CRAN](http://cran.r-project.org/web/packages/jaatha/index.html)