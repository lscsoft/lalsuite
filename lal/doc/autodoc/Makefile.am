## Process this file with automake to produce Makefile.in
MOSTLYCLEANFILES = logfile *.log *.ilg *.blg *.out
CLEANFILES = *.aux *.bbl *.dvi *.glo *.idx *.ind *.lof *.lot *.toc \
             *.eps *.pdf *.ps
DISTCLEANFILES = *.sty *.cls *.tex

dvi-local:
	rm -f *.aux
	@ docdir=`cd $(top_srcdir)/doc && pwd` ; \
          sources="$$docdir/ligodcc.cls $$docdir/lal.sty $$docdir/lalspec.tex $$docdir/lsd.tex $$docdir/lsdFig* $$docdir/schedule.*"; \
          for file in $$sources ; do \
            if test ! -f `basename $$file` ; then \
	      cmd="@LN_S@ $$file ."; echo "$$cmd"; $$cmd; \
	    fi; \
          done
	@LATEX@ lalspec && @LATEX@ lalspec && @LATEX@ lalspec
	pool_size=400000 @LATEX@ lsd && @MKIND@ lsd && \
	pool_size=400000 @LATEX@ lsd && \
	pool_size=400000 @LATEX@ lsd
	@ for file in lalspec.dvi lalspec.pdf lalspec.ps lsd.dvi lsd.pdf lsd.ps; do \
	    if test -f $$file ; then \
	      cmd="cp $$file .."; echo "$$cmd"; $$cmd; \
	    fi; \
	  done
