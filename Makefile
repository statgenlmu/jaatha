output = tmp
docPath = jaatha/vignettes
rnwfile = jaatha

default: package

howto:
	# Builds and opens Jaatha's vignette 
	- mkdir $(output)
	cd $(output); R CMD Sweave ../$(docPath)/$(rnwfile).Rnw
	- mv $(docPath)/*.tex $(docPath)/*.pdf $(docPath)/*.eps $(output)
	- cd $(docPath); cp *.png *.bib ../../$(output)/
	cd $(output); pdflatex $(rnwfile).tex
	cd $(output); bibtex $(rnwfile)
	cd $(output); pdflatex $(rnwfile).tex
	cd $(output); pdflatex $(rnwfile).tex
	cd $(output); evince $(rnwfile).pdf &

doc: clean-doc
	# Builds the roxygen2 documentation of Jaatha
	Rscript -e 'library(roxygen2); roxygenise("jaatha")'

test: doc
	# Runs the unit tests
	R CMD INSTALL jaatha
	cd unit_tests; ./doRUnit.R

check: doc clean-package
	# Runs an R CMD check
	R CMD check --no-vignettes jaatha
	make clean-package

package: test check
	# Build the R package out of the sources
	R CMD build jaatha

clean:
	# Deletes temoral output
	-rm -rv $(output) jaatha.Rcheck
	-rm -rv unitTests/results

clean-package:
	# Deletes temorary files, which tends to accumulate in the package
	- cd jaatha/src/; rm *.so *.o *.rds ms/*.o 2> /dev/null

clean-doc:
	# Deletes the automatically gernerated man files
	- rm jaatha/man/*.Rd 2> /dev/null
	- rm jaatha/DESCRIPTION 2> /dev/null 
	cp jaatha/DESCRIPTION.template jaatha/DESCRIPTION

commit: clean-package
	# Makes a git commit
	git add jaatha/R/*.R jaatha/src jaatha/inst/unitTests jaatha/DESCRIPTION Makefile
	git commit
	git push 
