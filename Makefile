.PHONY: howtos install test test-setup integration-test travis-test check clean

VERSION=$(shell grep Version DESCRIPTION | awk '{print $$2}')
PACKAGE=jaatha_$(VERSION).tar.gz

R_SOURCES=$(wildcard R/*.R) 
CPP_SOURCES=$(wildcard src/*.cc)
VIGNETTES=$(wildcard vignettes/*.pdf)
TESTS=$(wildcard inst/unitTests/*.R) $(wildcard tests/*.R)

default: $(PACKAGE)

release: clean test-setup howtos $(PACKAGE) check  
travis-test: $(PACKAGE) test-setup integration-test check

howtos: install 
	cd howtos; make
	cp howtos/*.pdf vignettes/

test: install
	cd tests; export RCMDCHECK=FALSE; Rscript doRUnit.R

integration-test: inst/unitTests/test_setup.Rda install 
	cd tests; export RCMDCHECK=FALSE; export INTEGRATION_TESTS=TRUE; Rscript doRUnit.R

test-setup: install
	cd inst/unitTests; Rscript test_setup.R

check: $(PACKAGE)
	# Runs an R CMD check
	R CMD check --as-cran $(PACKAGE)

package: $(PACKAGE) 

install: 
	R CMD INSTALL .

$(PACKAGE): $(R_SOURCES) $(CPP_SOURCES) $(TESTS) $(VIGNETTES) README DESCRIPTION man inst/unitTests/test_setup.Rda
	R CMD build .

README: README.md
	grep -v "\`\`\`" README.md | grep -v "Build Status" > README

inst/unitTests/test_setup.Rda: inst/unitTests/test_setup.R
	make test-setup

man: $(R_SOURCES) DESCRIPTION
	- rm -r man 2> /dev/null
	Rscript -e 'library(roxygen2); roxygenise(".")'

clean:
	- rm -rv jaatha.Rcheck
	- cd src/; rm *.so *.o *.rds ms/*.o 2> /dev/null
	- rm -r man 2> /dev/null
	- cd howtos; make clean
	- rm -rv inst/unitTests/test_setup.Rda
