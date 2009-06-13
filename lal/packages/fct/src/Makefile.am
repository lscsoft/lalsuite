# Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libfct.la
EXTRA_DIST = fct.h fct.c fct_fft.h fct_fft.c

if INTELFFT
libfct_la_SOURCES =
else
libfct_la_SOURCES = LALfct.c LALFCTInterface.c
endif

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
