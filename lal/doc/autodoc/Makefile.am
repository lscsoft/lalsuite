## Process this file with automake to produce Makefile.in
MOSTLYCLEANFILES = logfile *.log *.ilg *.blg *.out
CLEANFILES = *.aux *.bbl *.dvi *.glo *.idx *.ind *.lof *.lot *.tex *.toc \
             *.eps *.pdf *.ps

dvi-local:
	cp $(srcdir)/../lsd.* $(srcdir)
	rm -f lsdFig*
	@LN_S@ ../lsdFig* .
	rm -f schedule.*
	@LN_S@ ../schedule.* .
	@LATEX@ lsd.tex
	@MKIND@ lsd.idx
	@LATEX@ lsd.tex
	@LATEX@ lsd.tex
	rm -f $(srcdir)/lsd.tex
	mv $(srcdir)/lsd.* $(srcdir)/..
