all : sujet.pdf

sujet.pdf : sujet.tex                  \
            chapitres/introduction.tex \
            chapitres/objectifs.tex    \
            chapitres/preparation.tex  \
            chapitres/exercice.tex
	pdflatex sujet.tex ; pdflatex sujet.tex

show : sujet.pdf
	evince sujet.pdf&

clean : 
	rm -f *~ *.toc *.aux *.log *.lof *.bbl *.blg *.out 
	rm -f chapitres/*.aux                     chapitres/*~

mrproper : clean
	rm -f sujet.pdf
