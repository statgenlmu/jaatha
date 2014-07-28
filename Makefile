.PHONY: howtos install test test-setup integration-test travis-test check clean release

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

release: clean integration-test howtos $(PACKAGE) check

travis-test: $(PACKAGE) test-setup integration-test check

howtos: install 
	cd howtos; make

test: install
	export RCMDCHECK=FALSE; Rscript tests/doRUnit.R

integration-test: test-setup
	export RCMDCHECK=FALSE; export INTEGRATION_TESTS=TRUE; Rscript tests/doRUnit.R

test-setup: inst/unitTests/test-setup.Rda

check: $(PACKAGE)
	# check: Runs an R CMD check
	R CMD check $(R_CHECK_ARGS) $(PACKAGE)

package: $(PACKAGE) 
	# package: Builds the package

install: 
	# install: Installs the current version
	Rscript -e "library(Rcpp); compileAttributes()"
	Rscript -e 'library(roxygen2); roxygenise(".")'
	R CMD INSTALL .

clean:
	# clean: Removes temporary files
	- rm README NEWS NAMESPACE
	- rm -r jaatha.Rcheck 2> /dev/null
	- cd src/; rm *.so *.o *.rds ms/*.o 2> /dev/null
	- rm -r man 2> /dev/null
	- rm -r inst/unitTests/test-setup.Rda 2> /dev/null
	- cd howtos; make clean

#--------------------------------------------------------------------------
# Real targets
#--------------------------------------------------------------------------
inst/unitTests/test-setup.Rda: inst/unitTests/test-setup.R $(R_SOURCES) $(CPP_SOURCES)
	make install
	cd inst/unitTests; Rscript test-setup.R

$(PACKAGE): $(R_SOURCES) $(CPP_SOURCES) $(TESTS) README NEWS inst/unitTests/test-setup.Rda
	Rscript -e 'library(Rcpp); compileAttributes()'
	Rscript -e 'library(roxygen2); roxygenise(".")'
	R CMD build $(R_BUILD_ARGS) .

README: README.md
	grep -v "\`\`\`" README.md | grep -v "Build Status" > README

NEWS: NEWS.md
	cp NEWS.md NEWS
