## Process this file with automake to produce Makefile.in
MOSTLYCLEANFILES = logfile *.log *.ilg *.blg *.out
CLEANFILES = *.aux *.bbl *.dvi *.glo *.idx *.ind *.lof *.lot *.toc \
             *.eps *.pdf *.ps
DISTCLEANFILES = *.sty *.cls *.tex
TEXFLAGS = @TEXFLAGS@
LATEX = @LATEX@
MKIND = @MKIND@

dvi-local:
	rm -f *.aux
	sed '/^[A-Z]/s/$$/\\\\/' $(top_srcdir)/AUTHORS > AUTHORS.tex
	@ docdir="$(abs_top_srcdir)/doc" ; \
          sources="$$docdir/ligodcc.cls $$docdir/lal.sty $$docdir/lalspec.tex $$docdir/lsd.tex $$docdir/lsdFig* $$docdir/merger.*"; \
          for file in $$sources ; do \
            if test ! -f `basename $$file` ; then \
	      cmd="@LN_S@ $$file ."; echo "$$cmd"; $$cmd; \
	    fi; \
          done
	${LATEX} lalspec ${TEXFLAGS} && ${LATEX} lalspec ${TEXFLAGS} \
	  && ${LATEX} lalspec ${TEXFLAGS}
	max_strings=40000 pool_size=400000 ${LATEX} lsd ${TEXFLAGS} && \
	max_strings=40000 pool_size=400000 ${LATEX} lsd ${TEXFLAGS} && ${MKIND} lsd && \
	max_strings=40000 pool_size=400000 ${LATEX} lsd ${TEXFLAGS}
	@ for file in lalspec.dvi lalspec.pdf lalspec.ps lsd.dvi lsd.pdf lsd.ps; do \
	    if test -f $$file ; then \
	      cmd="cp $$file .."; echo "$$cmd"; $$cmd; \
	    fi; \
	  done
