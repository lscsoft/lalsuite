## Process this file with automake to produce Makefile.in
DOCSOURCEFILES = window.tex windowFig.eps windowFig.pdf \
  window_t.eps window_t.pdf window_f.eps window_f.pdf
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = main.dvi main.pdf main.ps
EXTRA_DIST = $(DOCSOURCEFILES)
dvi-local: .dvi-dep
distclean-local:
	-rm -fr .adoc
