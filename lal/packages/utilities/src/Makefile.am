## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libutilities.la
libutilities_la_SOURCES = \
  Random.c FindRoot.c Integrate.c Interpolate.c HeapSort.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
