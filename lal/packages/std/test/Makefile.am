## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALVersionTest LALMallocTest LALGSLTest StringConvertTest

LALVersionTest_SOURCES = LALVersionTest.c
LALVersionTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALMallocTest_SOURCES = LALMallocTest.c
LALMallocTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALGSLTest_SOURCES = LALGSLTest.c
LALGSLTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
StringConvertTest_SOURCES = StringConvertTest.c
StringConvertTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = LALVersionTest LALMallocTest LALGSLTest StringConvertTest
TESTS_ENVIRONMENT = 1>/dev/null

EXTRA_DIST = LALPrimerTest.c
EXTRA_TEXSOURCES = LALPrimerTest.c
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
