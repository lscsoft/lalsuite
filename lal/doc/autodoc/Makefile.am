## Process this file with automake to produce Makefile.in
MOSTLYCLEANFILES = logfile *.log *.ilg *.blg *.out
CLEANFILES = *.aux *.bbl *.dvi *.glo *.idx *.ind *.lof *.lot *.tex *.toc \
             *.eps *.ps

dvi-local:
	cp $(srcdir)/../lsd.* $(srcdir)
	@LATEX@ lsd.tex
	makeindex lsd.idx
	@LATEX@ lsd.tex
	@LATEX@ lsd.tex
	rm -f $(srcdir)/lsd.tex
	mv $(srcdir)/lsd.* $(srcdir)/..
