## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Window.h
CLEANFILES =
all-local: include-link
@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
CLEANFILES += $(swig_clean)
@SWIG_HEADER_MK@
