## Process this file with automake to produce Makefile.in
check_PROGRAMS = PrintVectorTest\
	         PrintFTSeriesTest ReadFTSeriesTest StreamInputTest\
                 StreamSeriesInputTest SFTfileIOTest \
		 ConfigFileTest UserInputTest LALMath3DPlotTest LALMathNDPlotTest \
		 SegmentsIOTest NRWaveIOTest LALXMLTest

noinst_DATA = StreamInput.data StreamSeriesInput.data inputsft.0 inputsft.1 ConfigFileSample.cfg SegmentsInput1.data \
	SFT-bad1   SFT-bad13  SFT-bad4	SFT-bad8   SFT-test2  SFT-test6 SFT-bad10  SFT-bad14  SFT-bad5	SFT-bad9   \
	SFT-test3  SFT-test7  SFT-bad11  SFT-bad2   SFT-bad6	SFT-good   SFT-test4  SFT-bad12  SFT-bad3   SFT-bad7	\
	SFT-test1  SFT-test5

EXTRA_DIST = $(DATA)
AM_CPPFLAGS = -DSRCDIR=\"$(srcdir)\"

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

LALMath3DPlotTest_SOURCES = LALMath3DPlotTest.c
LALMath3DPlotTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

LALMathNDPlotTest_SOURCES = LALMathNDPlotTest.c
LALMathNDPlotTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

SegmentsIOTest_SOURCES = SegmentsIOTest.c
SegmentsIOTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

NRWaveIOTest_SOURCES = NRWaveIOTest.c
NRWaveIOTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

LALXMLTest_SOURCES = LALXMLTest.c
LALXMLTest_LDADD = $(top_builddir)/packages/support/src/liblalxml.la $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
#TESTS = PrintVectorTest \
        #ReadFTSeriesTest PrintFTSeriesTest StreamInputTest \
        #StreamSeriesInputTest ConfigFileTest UserInputTest SFTfileIOTest \
        #LALMath3DPlotTest LALMathNDPlotTest
#TESTS = PrintVectorTest \
        #ReadFTSeriesTest PrintFTSeriesTest StreamInputTest \
        #StreamSeriesInputTest ConfigFileTest UserInputTest SFTfileIOTest \
        #LALMath3DPlotTest LALMathNDPlotTest SegmentsIOTest
TESTS = PrintVectorTest \
        ReadFTSeriesTest PrintFTSeriesTest StreamInputTest \
        StreamSeriesInputTest ConfigFileTest UserInputTest \
        LALMath3DPlotTest LALMathNDPlotTest SegmentsIOTest \
        LALXMLTest SFTfileIOTest

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null LAL_DATA_PATH=$(srcdir)
MOSTLYCLEANFILES = .dvi-dep *.dat *PrintVector.00* *.out  outputsft.0 \
		   TestOutputSFT.0 TestOutputSFT.1 Math3DNotebook.nb MathNDNotebook.nb \
		   SegmentsOutput1.data outputsft*.sft H-1_H1*.sft
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
