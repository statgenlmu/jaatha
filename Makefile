.PHONY: howtos install test quick-test check clean

VERSION=$(shell grep Version DESCRIPTION.template | awk '{print $$2}')
PACKAGE=jaatha_$(VERSION).tar.gz

R_SOURCES=$(wildcard R/*.R) 
CPP_SOURCES=$(wildcard src/*.cc)
TESTS=$(wildcard inst/unitTests/*.R)

default: $(PACKAGE)

release: $(PACKAGE) test check howtos 

howtos: install 
	cd howtos; make

test: install
	# Runs the unit tests
	cd tests; export RCMDCHECK=FALSE; ./doRUnit.R

quick-test: install
	# Runs the unit tests without time-consuming whole algorithms tests
	cd unit_tests; ./doRUnit.R quick

check: install 
	# Runs an R CMD check
	R CMD check --as-cran $(PACKAGE)

package: test check
	# Build the R package out of the sources
	R CMD build .

install: $(PACKAGE)
	R CMD INSTALL $(PACKAGE)

$(PACKAGE): $(R_SOURCES) $(CPP_SOURCES) $(TESTS) README DESCRIPTION man
	R CMD build .

README: README.md
	grep -v "\`\`\`" README.md > README

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
