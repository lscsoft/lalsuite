INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

check_PROGRAMS = \
	LALGSLTest \
	LALMallocTest \
	LALVersionTest \
	StringConvertTest

LALGSLTest_SOURCES = LALGSLTest.c
LALGSLTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

LALMallocTest_SOURCES = LALMallocTest.c
LALMallocTest_LDADD =  \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

LALVersionTest_SOURCES = LALVersionTest.c
LALVersionTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

StringConvertTest_SOURCES = StringConvertTest.c
StringConvertTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

TESTS = \
	LALGSLTest \
	LALMallocTest \
	LALVersionTest \
	StringConvertTest

TESTS_ENVIRONMENT = 1>/dev/null

EXTRA_DIST =

MOSTLYCLEANFILES =
DISTCLEANFILES =
