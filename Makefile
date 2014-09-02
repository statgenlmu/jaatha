.PHONY: howtos install test test-setup fulltest travis-test check clean release

VERSION=$(shell grep Version DESCRIPTION | awk '{print $$2}')
PACKAGE=jaatha_$(VERSION).tar.gz
R_CHECK_ARGS?="--as-cran"
R_BUILD_ARGS?=

R_SOURCES=$(wildcard R/*.R) 
CPP_SOURCES=$(wildcard src/*.cc)
TESTS=$(wildcard inst/unitTests/*.R) $(wildcard tests/*.R)

#--------------------------------------------------------------------------
# Phony targets
#--------------------------------------------------------------------------
default: $(PACKAGE)

release: clean fulltest howtos $(PACKAGE) check

travis-test: $(PACKAGE) test-setup fulltest check

howtos: install 
	cd howtos; make

test: install
	Rscript -e "library(devtools); test()"

integration-test: install 
	cd tests/integration; ./run-tests.R

check: $(PACKAGE)
	# check: Runs an R CMD check
	R CMD check $(R_CHECK_ARGS) $(PACKAGE)

package: $(PACKAGE) 
	# package: Builds the package

install: 
	# install: Installs the current version
	Rscript -e "library(Rcpp); compileAttributes()"
	Rscript -e 'library(roxygen2); roxygenise(".")'
	R CMD INSTALL --install-tests .

clean:
	# clean: Removes temporary files
	- rm README NEWS NAMESPACE
	- rm -r jaatha.Rcheck 2> /dev/null
	- cd src/; rm *.so *.o *.rds ms/*.o 2> /dev/null
	- rm -r man 2> /dev/null
	- cd howtos; make clean

#--------------------------------------------------------------------------
# Real targets
#--------------------------------------------------------------------------

$(PACKAGE): $(R_SOURCES) $(CPP_SOURCES) $(TESTS) README NEWS
	Rscript -e 'library(Rcpp); compileAttributes()'
	Rscript -e 'library(roxygen2); roxygenise(".")'
	R CMD build $(R_BUILD_ARGS) .

README: README.md
	grep -v "\`\`\`" README.md | grep -v "Build Status" > README

NEWS: NEWS.md
	cp NEWS.md NEWS
