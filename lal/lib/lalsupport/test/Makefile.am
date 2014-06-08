include $(top_srcdir)/gnuscripts/lallibs_header_links.am

AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include -DTEST_DATA_DIR='"$(abs_srcdir)/"' @AM_CPPFLAGS@
LDADD = \
	$(top_builddir)/lib/lal/liblal.la \
	$(top_builddir)/lib/lalsupport/src/liblalsupport.la

check_PROGRAMS = \
	ConfigFileTest \
	LALMath3DPlotTest \
	LALMathNDPlotTest \
	PrintFTSeriesTest \
	PrintVectorTest \
	ReadFTSeriesTest \
	SegmentsIOTest \
	StreamInputTest \
	StreamSeriesInputTest \
	UserInputTest \
	GzipTest

noinst_DATA = \
	ConfigFileSample.cfg \
	ConfigFileSample2.cfg \
	ConfigFileSample3.cfg.gz \
	SegmentsInput1.data \
	StreamInput.data \
	StreamSeriesInput.data

EXTRA_DIST = \
	$(DATA) \
	GzipTest.sh

TESTS = \
	$(check_PROGRAMS) \
	GzipTest.sh

MOSTLYCLEANFILES = \
	*.dat \
	*.out \
	*PrintVector.00* \
	Math3DNotebook.nb \
	MathNDNotebook.nb \
	SegmentsOutput1.data \
	ConfigFile.cfg
