AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include -DPKG_DATA_DIR='"$(pkgdatadir)/"' @AM_CPPFLAGS@

lib_LTLIBRARIES = liblalsupport.la

liblalsupport_la_SOURCES = \
	Audio.c \
	ConfigFile.c \
	FileIO.c \
	LALCache.c \
	LALMath3DPlot.c \
	LALMathNDPlot.c \
	LogPrintf.c \
	PrintFrequencySeries.c \
	PrintTimeSeries.c \
	PrintVector.c \
	ReadFiltersFile.c \
	ReadFrequencySeries.c \
	ReadTimeSeries.c \
	ReadNoiseSpectrum.c \
	SegmentsIO.c \
	StreamSeriesInput.c \
	StreamSeriesOutput.c \
	StreamSequenceInput.c \
	StreamGridInput.c \
	StreamGridOutput.c \
	StreamVectorInput.c \
	StreamVectorSequenceInput.c \
	UserInput.c \
	LALgetopt.c \
	LALgetopt.h


liblalsupport_la_LIBADD = $(top_builddir)/lib/lal/liblal.la $(ZLIB_LIBS)
liblalsupport_la_CFLAGS = $(AM_CFLAGS) $(ZLIB_CFLAGS)
liblalsupport_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION_SUPPORT)

noinst_HEADERS = \
	PrintFrequencySeries_source.c \
	PrintTimeSeries_source.c \
	PrintVector_source.c \
	ReadFrequencySeries_source.c \
	ReadTimeSeries_source.c \
	StreamGridInput_source.c \
	StreamGridOutput_source.c \
	StreamSequenceInput_source.c \
	StreamSequenceInputComplex_source.c \
	StreamSeriesInput_source.c \
	StreamSeriesOutput_source.c \
	StreamVectorInput_source.c \
	StreamVectorSequenceInput_source.c
