## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALMomentTest

LALMomentTest_SOURCES = LALMomentTest.c CheckStatus.c
LALMomentTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = LALMomentTest

EXTRA_DIST = CheckStatus.c CheckStatus.h 

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
