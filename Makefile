output = tmp
docPath = package/vignettes
rnwfile = jaatha

default: package

howtos: 
	# Builds and opens Jaatha's vignette 
	- mkdir $(output) 2> /dev/null;
	- mkdir $(docPath)/cache 2> /dev/null;
	cp -r $(docPath)/* $(output)/
	cd $(output);\
	  ../misc/knitr.R ../$(docPath)/custom_simulator_howto.Rnw ../$(docPath)/cache/
	#cd $(output); ../misc/knitr.R ../$(docPath)/jaatha_howto.Rnw ../$(docPath)/cache/

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

package: test check 
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

package/vignettes/jaatha_ascii.bib: package/vignettes/jaatha_utf8.bib
	cd package/vignettes; iconv -f=utf8 -t=ascii jaatha.utf8.bib > jaatha_ascii.bib
