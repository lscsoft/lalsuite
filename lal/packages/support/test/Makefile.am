## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALVersionTest LALMallocTest LALGSLTest PrintVectorTest\
	         PrintFTSeriesTest ReadFTSeriesTest StreamInputTest\
                 StreamSeriesInputTest SFTfileIOTest \
		ConfigFileTest UserInputTest

noinst_DATA = StreamInput.data StreamSeriesInput.data inputsft.0
EXTRA_DIST = $(DATA)
LALVersionTest_SOURCES = LALVersionTest.c
LALVersionTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALMallocTest_SOURCES = LALMallocTest.c
LALMallocTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALGSLTest_SOURCES = LALGSLTest.c
LALGSLTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

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

SFTfileIOTest_SOURCES = SFTfileIOTest.c
SFTfileIOTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

ConfigFileTest_SOURCES = ConfigFileTest.c
ConfigFileTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

UserInputTest_SOURCES = UserInputTest.c
UserInputTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = LALVersionTest LALMallocTest LALGSLTest PrintVectorTest \
        ReadFTSeriesTest PrintFTSeriesTest StreamInputTest \
        StreamSeriesInputTest ConfigFileTest UserInputTest SFTfileIOTest

TESTS_ENVIRONMENT = LAL_DATA_PATH=$(srcdir)
MOSTLYCLEANFILES = .dvi-dep *.dat *PrintVector.00* *.out  outputsft.0 
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
