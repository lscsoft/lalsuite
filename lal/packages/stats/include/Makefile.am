## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALMoment.h LALCorrelation.h
CLEANFILES =
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
all-local: include-link
dvi-local: .dvi-dep
@INCLUDEMKFILE@
@DVIMKFILE@

swig_headers = $(pkginclude_HEADERS)
CLEANFILES += $(swig_clean)
@SWIG_HEADER_MK@
