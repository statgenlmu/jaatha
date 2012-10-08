## Makefile to include unit tests on R CMD check
## Taken from:
## http://rwiki.sciviews.org/doku.php?id=developers:runit
## See there for copyright terms

TOP=".."
PKG="jaatha"

all: inst test

inst: # Install package
	cd ${TOP};\
	R CMD INSTALL ${PKG}

test: # Run unit tests
	./doRUnit.R
