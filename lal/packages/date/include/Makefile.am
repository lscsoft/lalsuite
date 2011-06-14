## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Date.h TimeDelay.h
if SWIG_BUILD
pkginclude_HEADERS += swig-lal-date.i
endif
CLEANFILES =
MOSTLYCLEANFILES =
DISTCLEANFILES =
all-local: include-link
@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
CLEANFILES += $(swig_clean)
@SWIG_HEADER_MK@
