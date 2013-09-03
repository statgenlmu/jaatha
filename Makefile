output = tmp
docPath = package/vignettes
rnwfile = jaatha

default: package

howto:
	# Builds and opens Jaatha's vignette 
	- mkdir $(output); cp -r $(docPath)/* $(output)/
	cd $(output); R CMD Sweave $(rnwfile).Rnw;
	cd $(output); pdflatex $(rnwfile).tex;
	cd $(output); bibtex $(rnwfile);
	cd $(output); pdflatex $(rnwfile).tex;
	cd $(output); pdflatex $(rnwfile).tex;
	cd $(output); evince $(rnwfile).pdf &

howto-cache:
	- cd $(docPath); R CMD Sweave jaatha.Rnw
	cd $(docPath)/cache; R CMD Sweave initialSearch.Rnw
	#- cd $(docPath); R CMD Sweave jaatha.Rnw
	cd $(docPath)/cache; R CMD Sweave refineSearch.Rnw
	rm $(docPath)/jaatha.tex

doc: clean-doc
	# Builds the roxygen2 documentation of Jaatha
	Rscript -e 'library(roxygen2); roxygenise("package")'

test: doc
	# Runs the unit tests
	R CMD INSTALL package
	cd unit_tests; ./doRUnit.R

quick-test: doc
	# Runs the unit tests without time-consuming whole algorithms tests
	R CMD INSTALL package
	cd unit_tests; ./doRUnit.R quick

check: doc clean-package
	# Runs an R CMD check
	R CMD check package
	make clean-package

package: test howto-cache check
	# Build the R package out of the sources
	R CMD build package

clean:
	# Deletes temoral output
	-rm -rv $(output) jaatha.Rcheck
	-rm -rv unitTests/results

clean-package:
	# Deletes temorary files, which tends to accumulate in the package
	- cd package/src/; rm *.so *.o *.rds ms/*.o 2> /dev/null

clean-doc:
	# Deletes the automatically gernerated man files
	- rm package/man/*.Rd 2> /dev/null
	- rm package/DESCRIPTION 2> /dev/null 
	cp DESCRIPTION.template package/DESCRIPTION
