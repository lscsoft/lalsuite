## Process this file with automake to produce Makefile.in
DOCSOURCEFILES = tools.tex \
  resamp_1000x4.eps resamp_1000x4.pdf resamp_1000x8.eps resamp_1000x8.pdf \
  resamp_100x2_end.eps resamp_100x2_end.pdf resamp_100x2_start.eps \
  resamp_100x2_start.pdf resamp_100x4.eps resamp_100x4.pdf resamp_10x4.eps \
  resamp_10x4.pdf \
  sequence.tex frequencyseries.tex timeseries.tex
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = main.dvi main.pdf main.ps
EXTRA_DIST = $(DOCSOURCEFILES)
dvi-local: .dvi-dep
distclean-local:
	-rm -fr .adoc
@DOCMKFILE@
