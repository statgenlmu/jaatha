## Makefile to include unit tests on R CMD check
## Taken from:
## http://rwiki.sciviews.org/doku.php?id=developers:runit
## See there for copyright terms

TOP=../..
PKG=${shell cd ${TOP};pwd}
SUITE=doRUnit.R
R=R

all: inst test

inst: # Install package
	cd ${TOP}/..;\
	${R} CMD INSTALL ${PKG}

test: # Run unit tests
	export RCMDCHECK=FALSE;\
	cd ${TOP}/tests;\
	${R} --vanilla --slave < ${SUITE}
