## Process this file with automake to produce Makefile.in
MOSTLYCLEANFILES = logfile *.log *.ilg *.blg *.out
CLEANFILES = *.aux *.bbl *.dvi *.glo *.idx *.ind *.lof *.lot *.tex *.toc \
             *.eps *.pdf *.ps

dvi-local:
	cp $(srcdir)/../ligodcc.cls .
	cp $(srcdir)/../lalspec.* .
	@LATEX@ lalspec.tex
	@LATEX@ lalspec.tex
	@LATEX@ lalspec.tex
	mv lalspec.* ..
	rm -f ligodcc.cls
	cp $(srcdir)/../lsd.* .
	rm -f lsdFig*
	@LN_S@ `pwd`/$(srcdir)/../lsdFig* .
	rm -f schedule.*
	@LN_S@ `pwd`/$(srcdir)/../schedule.* .
	@LATEX@ lsd.tex
	@MKIND@ lsd.idx
	@LATEX@ lsd.tex
	@LATEX@ lsd.tex
	rm -f lsd.tex
	mv lsd.* ..
