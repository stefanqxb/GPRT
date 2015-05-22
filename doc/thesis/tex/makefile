all: xjobb-en.pdf


%.pdf: %.tex
	pdflatex $*.tex
	if ( grep -q citation $*.aux ) ; then \
		bibtex $* ; \
		pdflatex $*.tex ; \
	fi;
	pdflatex $*.tex

rem:
	rm thesis.pdf
clean:    
	find . -regex ".*\.\(aux\|bbl\|blg\|log\|lof\|toc\|lot\|nlo\|out\)" -exec rm {} \;

