## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libstats.la
libstats_la_SOURCES = LALMoment.c LALCorrelation.c
noinst_HEADERS = LALMoment_source.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

@DVIMKFILE@
