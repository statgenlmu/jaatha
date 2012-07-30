output = tmp
docPath = jaatha/vignettes
rnwfile = jaatha

default: package

howto:
	- mkdir $(output)
	cd $(output); R CMD Sweave ../$(docPath)/$(rnwfile).Rnw
	- mv $(docPath)/*.tex $(docPath)/*.pdf $(docPath)/*.eps $(output)
	- cd $(docPath); cp *.png *.bib ../../$(output)/
	cd $(output); pdflatex $(rnwfile).tex
	cd $(output); bibtex $(rnwfile)
	cd $(output); pdflatex $(rnwfile).tex
	cd $(output); pdflatex $(rnwfile).tex
	cd $(output); evince $(rnwfile).pdf &

doc:
	Rscript -e 'library(roxygen2); roxygenise("jaatha")'

test:
	cd jaatha/inst/unitTests; make

check: doc
	R CMD check --no-vignettes jaatha
	make clean-package

package: check
	R CMD build jaatha

clean:
	-rm -rv $(output) jaatha.Rcheck

clean-package:
	- cd jaatha/inst/unitTests; rm -v report* 2> /dev/null
	- cd jaatha/src/; rm -v *.so *.o *.rds ms/*.o 2> /dev/null

commit: clean-package
	git add jaatha/R/*.R jaatha/src jaatha/inst/unitTests jaatha/DESCRIPTION Makefile
	git commit
	git push neon master
