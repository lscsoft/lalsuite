## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libring.la
libring_la_SOURCES = Ring.c RingSearchInit.c RingSearchConditionData.c \
    RingSearch.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
