## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALVersionTest LALMallocTest LALGSLTest

LALVersionTest_SOURCES = LALVersionTest.c
LALVersionTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALMallocTest_SOURCES = LALMallocTest.c
LALMallocTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALGSLTest_SOURCES = LALGSLTest.c
LALGSLTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = LALVersionTest LALMallocTest LALGSLTest

EXTRA_DIST = LALPrimerTest.c
EXTRA_TEXSOURCES = LALPrimerTest.c
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
