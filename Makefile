all: decomp.pdf

decomp.pdf: decomp.tex
	pdflatex decomp.tex

decomp.tex: decomp.raw.tex
	./insert_code.pl decomp.g <decomp.raw.tex >decomp.tex

.PHONY: clean
clean:
	rm -f decomp.tex *.pdf *.thm *.log *.toc
