output = tmp
docPath = jaatha/vignettes
rnwfile = jaatha

doc:
	    - mkdir $(output)
	     cd $(output); R CMD Sweave ../$(docPath)/$(rnwfile).Rnw
	    - mv $(docPath)/*.tex $(docPath)/*.pdf $(docPath)/*.eps $(output)
	    - cd $(docPath); cp *.png ../../$(output)/
	    cd $(output); pdflatex $(rnwfile).tex
	    cd $(output); pdflatex $(rnwfile).tex
	    cd $(output); evince $(rnwfile).pdf &

clean:
	    -rm $(output)
