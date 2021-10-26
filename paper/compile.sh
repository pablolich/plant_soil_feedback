pdflatex paper.tex
pdflatex paper.tex
bibtex paper
pdflatex paper.tex
pdflatex paper.tex
evince paper.pdf &

pdflatex supp_methods.tex
pdflatex supp_methods.tex
bibtex supp_methods
pdflatex supp_methods.tex
pdflatex supp_methods.tex
evince supp_methods.pdf &

rm *.aux *.bbl *.blg *.log *.out *.synctex.gz
