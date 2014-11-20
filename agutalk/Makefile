all: talk.pdf

talk.aux: talk.tex
	pdflatex talk

talk.pdf: talk.aux
	pdflatex talk

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *.nav *.snm *.toc *~
