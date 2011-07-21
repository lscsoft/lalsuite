## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALMoment.h LALCorrelation.h
CLEANFILES =
MOSTLYCLEANFILES =
DISTCLEANFILES =
all-local: include-link

@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
CLEANFILES += $(swig_clean)
@SWIG_HEADER_MK@
