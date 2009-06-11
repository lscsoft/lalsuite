## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libclremoval.la
libclremoval_la_SOURCES = CleanAll.c HarmonicFinder.c RefInterference.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
