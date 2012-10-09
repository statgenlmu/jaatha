output = tmp
docPath = jaatha/vignettes
rnwfile = jaatha

default: package

howto:
	# Builds and opens Jaatha's vignette 
	- mkdir $(output); cp -r $(docPath)/* $(output)/
	cd $(output); R CMD Sweave $(rnwfile).Rnw;\
				  pdflatex $(rnwfile).tex;\
				  bibtex $(rnwfile);\
				  pdflatex $(rnwfile).tex;\
				  pdflatex $(rnwfile).tex;\
				  evince $(rnwfile).pdf &

howto-cache:
	cd $(docPath); R CMD Sweave jaatha.Rnw
	cd $(docPath)/cache; R CMD Sweave initialSearch.Rnw
	cd $(docPath); R CMD Sweave jaatha.Rnw
	cd $(docPath)/cache; R CMD Sweave refineSearch.Rnw
	rm $(docPath)/jaatha.tex

doc: clean-doc
	# Builds the roxygen2 documentation of Jaatha
	Rscript -e 'library(roxygen2); roxygenise("jaatha")'

test: doc
	# Runs the unit tests
	cd unit_tests; ./doRUnit.R

check: doc clean-package
	# Runs an R CMD check
	R CMD check jaatha
	make clean-package

package: test howto-cache check
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
