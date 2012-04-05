INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

check_PROGRAMS = \
	ConfigFileTest \
	ConfigFileTest2 \
	LALMath3DPlotTest \
	LALMathNDPlotTest \
	PrintFTSeriesTest \
	PrintVectorTest \
	ReadFTSeriesTest \
	SegmentsIOTest \
	StreamInputTest \
	StreamSeriesInputTest \
	UserInputTest

noinst_DATA = \
	ConfigFileSample.cfg \
	ConfigFileSample2.cfg \
	SegmentsInput1.data \
	StreamInput.data \
	StreamSeriesInput.data

EXTRA_DIST = $(DATA)

ConfigFileTest_SOURCES = ConfigFileTest.c
ConfigFileTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

ConfigFileTest2_SOURCES = ConfigFileTest2.c
ConfigFileTest2_LDADD = \
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
	ConfigFileTest2 \
	LALMath3DPlotTest \
	LALMathNDPlotTest \
	PrintFTSeriesTest \
	PrintVectorTest \
	ReadFTSeriesTest \
	SegmentsIOTest \
	StreamInputTest \
	StreamSeriesInputTest \
	UserInputTest

MOSTLYCLEANFILES = \
	*.dat \
	*.out \
	*PrintVector.00* \
	Math3DNotebook.nb \
	MathNDNotebook.nb \
	SegmentsOutput1.data \
	ConfigFile.cfg

DISTCLEANFILES =
