AM_CPPFLAGS = -DSRCDIR=\"$(srcdir)\"
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

check_PROGRAMS = \
	ConfigFileTest \
	LALMath3DPlotTest \
	LALMathNDPlotTest \
	NRWaveIOTest \
	PrintFTSeriesTest \
	PrintVectorTest \
	ReadFTSeriesTest \
	SFTfileIOTest \
	SegmentsIOTest \
	StreamInputTest \
	StreamSeriesInputTest \
	UserInputTest

noinst_DATA = \
	ConfigFileSample.cfg \
	SFT-bad1 \
	SFT-bad10 \
	SFT-bad11 \
	SFT-bad12 \
	SFT-bad13 \
	SFT-bad14 \
	SFT-bad2 \
	SFT-bad3 \
	SFT-bad4 \
	SFT-bad5 \
	SFT-bad6 \
	SFT-bad7 \
	SFT-bad8 \
	SFT-bad9 \
	SFT-good \
	SFT-test1 \
	SFT-test2 \
	SFT-test3 \
	SFT-test4 \
	SFT-test5 \
	SFT-test6 \
	SFT-test7 \
	SegmentsInput1.data \
	StreamInput.data \
	StreamSeriesInput.data \
	inputsft.0 \
	inputsft.1

EXTRA_DIST = $(DATA)

ConfigFileTest_SOURCES = ConfigFileTest.c
ConfigFileTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

LALMath3DPlotTest_SOURCES = LALMath3DPlotTest.c
LALMath3DPlotTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

LALMathNDPlotTest_SOURCES = LALMathNDPlotTest.c
LALMathNDPlotTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

NRWaveIOTest_SOURCES = NRWaveIOTest.c
NRWaveIOTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

PrintFTSeriesTest_SOURCES = PrintFTSeriesTest.c
PrintFTSeriesTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

PrintVectorTest_SOURCES = PrintVectorTest.c
PrintVectorTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

ReadFTSeriesTest_SOURCES = ReadFTSeriesTest.c
ReadFTSeriesTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

SFTfileIOTest_SOURCES = SFTfileIOTest.c
SFTfileIOTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

SegmentsIOTest_SOURCES = SegmentsIOTest.c
SegmentsIOTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

StreamInputTest_SOURCES = StreamInputTest.c
StreamInputTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

StreamSeriesInputTest_SOURCES = StreamSeriesInputTest.c
StreamSeriesInputTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

UserInputTest_SOURCES = UserInputTest.c
UserInputTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null LAL_DATA_PATH=$(srcdir)

TESTS = \
	ConfigFileTest \
	LALMath3DPlotTest \
	LALMathNDPlotTest \
	PrintFTSeriesTest \
	PrintVectorTest \
	ReadFTSeriesTest \
	SFTfileIOTest \
	SegmentsIOTest \
	StreamInputTest \
	StreamSeriesInputTest \
	UserInputTest

MOSTLYCLEANFILES = \
	*.dat \
	*.out \
	*PrintVector.00* \
	.dvi-dep \
	H-1_H1*.sft \
	Math3DNotebook.nb \
	MathNDNotebook.nb \
	SegmentsOutput1.data \
	TestOutputSFT.0 \
	TestOutputSFT.1 \
	outputsft*.sft \
	outputsft.0

DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

@DVIMKFILE@
