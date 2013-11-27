.PHONY: howtos install test quick-test travis-test check clean

VERSION=$(shell grep Version DESCRIPTION.template | awk '{print $$2}')
PACKAGE=jaatha_$(VERSION).tar.gz
R_CHECK_ARGS?="--as-cran"
R_BUILD_ARGS?=""

R_SOURCES=$(wildcard R/*.R) 
CPP_SOURCES=$(wildcard src/*.cc)

default: $(PACKAGE)

release: $(PACKAGE) test check howtos 
travis-test: $(PACKAGE) test check

howtos: install 
	cd howtos; make

test: install
	# Runs the unit tests
	cd unit_tests; ./doRUnit.R

quick-test: install
	# Runs the unit tests without time-consuming whole algorithms tests
	cd unit_tests; ./doRUnit.R quick

check: install 
	# Runs an R CMD check
	R CMD check $(R_CHECK_ARGS) $(PACKAGE)

package: test check
	# Build the R package out of the sources
	R CMD build $(R_BUILD_ARGS) .

install:
	R CMD INSTALL .

$(PACKAGE): $(R_SOURCES) $(CPP_SOURCES) README DESCRIPTION man
	R CMD build .


README: README.md
	grep -v "\`\`\`" README.md | grep -v "Build Status" > README

DESCRIPTION: DESCRIPTION.template 
	cp DESCRIPTION.template DESCRIPTION

man: $(R_SOURCES) DESCRIPTION
	- rm -r man 2> /dev/null
	Rscript -e 'library(roxygen2); roxygenise(".")'

clean:
	- rm -rv jaatha.Rcheck
	- rm -rv unit_tests/results
	- cd src/; rm *.so *.o *.rds ms/*.o 2> /dev/null
	- rm -v DESCRIPTION README 2>/dev/null
	- rm -r man 2> /dev/null
