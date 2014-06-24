.PHONY: howtos install test test-setup integration-test travis-test check clean

VERSION=$(shell grep Version DESCRIPTION | awk '{print $$2}')
PACKAGE=jaatha_$(VERSION).tar.gz
R_CHECK_ARGS?="--as-cran"
R_BUILD_ARGS?=

R_SOURCES=$(wildcard R/*.R) 
CPP_SOURCES=$(wildcard src/*.cc)
VIGNETTES=$(wildcard vignettes/*.pdf)
TESTS=$(wildcard inst/unitTests/*.R) $(wildcard tests/*.R)

default: $(PACKAGE)

release: clean howtos $(PACKAGE) check  
travis-test: $(PACKAGE) test-setup integration-test check

howtos: install 
	cd howtos; make

test: install
	cd tests; export RCMDCHECK=FALSE; Rscript doRUnit.R

integration-test: inst/unitTests/test-setup.Rda install 
	cd tests; export RCMDCHECK=FALSE; export INTEGRATION_TESTS=TRUE; Rscript doRUnit.R

test-setup: install
	cd inst/unitTests; Rscript test-setup.R

check: $(PACKAGE)
	# Runs an R CMD check
	R CMD check $(R_CHECK_ARGS) $(PACKAGE)

package: $(PACKAGE) 

install: 
	Rscript -e "library(Rcpp); compileAttributes()"
	R CMD INSTALL .

$(PACKAGE): $(R_SOURCES) $(CPP_SOURCES) $(TESTS) $(VIGNETTES) README NEWS man test-setup
	R CMD build $(R_BUILD_ARGS) .

README: README.md
	grep -v "\`\`\`" README.md | grep -v "Build Status" > README

NEWS: NEWS.md
	cp NEWS.md NEWS

inst/unitTests/test-setup.Rda: inst/unitTests/test-setup.R
	make test-setup

man: $(R_SOURCES) DESCRIPTION
	- rm -r man 2> /dev/null
	Rscript -e 'library(roxygen2); roxygenise(".")'

clean:
	- rm README NEWS
	- rm -r jaatha.Rcheck 2> /dev/null
	- cd src/; rm *.so *.o *.rds ms/*.o 2> /dev/null
	- rm -r man 2> /dev/null
	- cd howtos; make clean
	- rm -r inst/unitTests/test-setup.Rda 2> /dev/null
