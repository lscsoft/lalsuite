## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALMomentTest LALCorrelationTest

LALMomentTest_SOURCES = LALMomentTest.c CheckStatus.c
LALMomentTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

LALCorrelationTest_SOURCES = LALCorrelationTest.c
LALCorrelationTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = LALMomentTest LALCorrelationTest
TESTS_ENVIRONMENT = 1>/dev/null

EXTRA_DIST = CheckStatus.c CheckStatus.h 

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
