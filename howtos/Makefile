howtos = jaatha_howto.pdf custom_simulator_howto.pdf

all: $(howtos)

%.tex: %.Rnw
	Rscript -e "library('knitr'); knit('$<')"; 

%.pdf: %.tex jaatha.bib
	pdflatex $<
	bibtex $(<:%.tex=%) 
	pdflatex $<
	pdflatex $<

clean:
	rm -r cache *.aux *.bbl *.blg *.log *.out *.tex *.pdf 2> /dev/null; true
