AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include
LDADD = \
	$(top_builddir)/lib/lal/liblal.la \
	$(top_builddir)/lib/lalsupport/src/liblalsupport.la

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

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null LAL_DATA_PATH=$(srcdir)

TESTS = \
	$(check_PROGRAMS)

MOSTLYCLEANFILES = \
	*.dat \
	*.out \
	*PrintVector.00* \
	Math3DNotebook.nb \
	MathNDNotebook.nb \
	SegmentsOutput1.data \
	ConfigFile.cfg
