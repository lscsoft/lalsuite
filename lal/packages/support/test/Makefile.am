## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALVersionTest LALMallocTest PrintVectorTest\
	         PrintFTSeriesTest ReadFTSeriesTest StreamInputTest\
                 StreamSeriesInputTest
noinst_DATA = StreamInput.data
EXTRA_DIST = $(DATA)
LALVersionTest_SOURCES = LALVersionTest.c
LALVersionTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALMallocTest_SOURCES = LALMallocTest.c
LALMallocTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

PrintVectorTest_SOURCES = PrintVectorTest.c
PrintVectorTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

PrintFTSeriesTest_SOURCES = PrintFTSeriesTest.c
PrintFTSeriesTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

ReadFTSeriesTest_SOURCES = ReadFTSeriesTest.c
ReadFTSeriesTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

StreamInputTest_SOURCES = StreamInputTest.c
StreamInputTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

StreamSeriesInputTest_SOURCES = StreamSeriesInputTest.c
StreamSeriesInputTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = LALVersionTest LALMallocTest PrintVectorTest \
        ReadFTSeriesTest PrintFTSeriesTest StreamInputTest \
        StreamSeriesInputTest
TESTS_ENVIRONMENT = LAL_DATA_PATH=$(srcdir)
MOSTLYCLEANFILES = .dvi-dep *.dat *PrintVector.00* *.out
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
