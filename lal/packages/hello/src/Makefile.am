## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libhello.la
libhello_la_SOURCES = LALHello.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
