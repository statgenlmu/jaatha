jaatha
======

Jaatha is a fast parameter estimation method primarily designed for Evolutionary
Biology. Details are describe in

> Mathew, L. A., Staab, P. R., Rose, L. E. and Metzler, D. (2013), 
[Why to account for finite sites in population genetic studies and how to do this with Jaatha 2.0](http://onlinelibrary.wiley.com/doi/10.1002/ece3.722/abstract). 
Ecology and Evolution.

Instructions for a quick-start are provided in the 
[The Jaatha HowTo](http://evol.bio.lmu.de/_statgen/software/jaatha/jaatha_howto.pdf). 
Instructions how to use Jaatha with a non-standard simulation method are given
in the 
[Custom Simulation Method HowTo](http://evol.bio.lmu.de/_statgen/software/jaatha/custom_simulator_howto.pdf).

# Installation
## From CRAN:
To install the current stable version of jaatha from cran, simply type

```R
install.packages('jaatha')
```

in R.

## Development version from GitHub  
You can download and build the development version from GitHub following the
commands: 

```bash
# Install Jaathas dependencies from CRAN
Rscript -e "install.packages(c('Rcpp', 'ape', 'doMC', 'RUnit'))"

# Download Jaatha from GitHub
git clone https://github.com/paulstaab/jaatha.git
cd jaatha

# And build the package using the Makefile
make
```

Afterwards you can install the created package with 
```bash
R CMD INSTALL jaatha_x.y.tar.gz
```

# Usage
Please refer to the 
[The Jaatha HowTo](http://evol.bio.lmu.de/_statgen/software/jaatha/jaatha_howto.pdf)
for usage information.

# Links
* (Jaatha's Homepage)[http://evol.bio.lmu.de/_statgen/software/jaatha]
* (Source Code on GitHub)[https://github.com/paulstaab/jaatha]
* (Bug tracker)[https://github.com/paulstaab/jaatha/issues]
* (Jaatha's page on CRAN)[http://cran.r-project.org/web/packages/jaatha/index.html]
