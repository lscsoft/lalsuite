## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libring.la
libring_la_SOURCES = Ring.c RingSearchInit.c RingSearchConditionData.c \
    RingSearch.c LALRingDown.c LALWiener.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
EXTRA_DIST = LALRingDown.m4 LALRingDownBase.m4

LALRingDown.c: LALRingDown.m4 LALRingDownBase.m4
	m4 LALRingDown.m4 > LALRingDown.c
