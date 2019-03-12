filename=prestocolor

pdf: pdf
	pdflatex ${filename}.tex; bibtex ${filename}; pdflatex ${filename}.tex; pdflatex ${filename}.tex


html:
	@#latex2html -split +0 -info "" -no_navigation ${filename}
	htlatex ${filename}

ps:	dvi
	dvips -t letter ${filename}.dvi

dvi:
	latex ${filename}
	bibtex ${filename}||true
	latex ${filename}
	latex ${filename}

clean:
	rm -f ${filename}.{log,aux,out,dvi,bbl,blg,4ct,4tc,css,idv,lg,tmp,xref,png,0x.png}

nuke:
	make clean
	rm -f ${filename}.{pdf,html}

