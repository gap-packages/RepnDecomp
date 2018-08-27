all: decomp.pdf

decomp.pdf: decomp.tex
	pdflatex decomp.tex
	pdflatex decomp.tex

decomp.tex: decomp.raw.tex
	./insert_code.pl ../gap/RepnDecomp.gi <decomp.raw.tex >decomp.tex

.PHONY: clean
clean:
	rm -f decomp.tex *.pdf *.thm *.log *.toc
