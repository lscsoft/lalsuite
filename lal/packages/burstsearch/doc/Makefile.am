## Process this file with automake to produce Makefile.in
DOCSOURCEFILES = burstsearch.tex
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = main.dvi main.pdf main.ps
EXTRA_DIST = $(DOCSOURCEFILES)
dvi-local: .dvi-dep
distclean-local:
	-rm -fr .adoc
@DOCMKFILE@
