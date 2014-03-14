all : nbody.pdf

nbody.pdf : nbody.tex                  \
            chapitres/introduction.tex \
            chapitres/objectifs.tex    \
            chapitres/preparation.tex  \
            chapitres/exercices.tex
	pdflatex nbody.tex ; pdflatex nbody.tex

show : nbody.pdf
	evince nbody.pdf&

clean : 
	rm -f *~ *.toc *.aux *.log *.lof *.bbl *.blg *.out 
	rm -f chapitres/*.aux                     chapitres/*~

mrproper : clean
	rm -f nbody.pdf
