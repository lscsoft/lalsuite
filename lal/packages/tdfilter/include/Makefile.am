## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LPC.h BandPassTimeSeries.h IIRFilter.h ZPGFilter.h
CLEANFILES =
MOSTLYCLEANFILES =
DISTCLEANFILES =
all-local: include-link
@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
