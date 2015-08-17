FILE = rtgp

pdf: *.tex
	pdflatex -shell-escape $(FILE)
bib:
	bibtex $(FILE)
clean: 
	-rm ./*.aux ./*.dvi ./*.blg ./*.log ./*.out ./*.toc ./*.spl ./*.tpt ./*.bbl ./$(FILE).pdf ./*~
