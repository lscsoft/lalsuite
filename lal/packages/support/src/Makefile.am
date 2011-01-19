lib_LTLIBRARIES = liblalsupport.la
liblalsupport_la_SOURCES = \
	Audio.c \
	ConfigFile.c \
	FileIO.c \
	LALCache.c \
	LALInitBarycenter.c \
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
	getopt.c \
	getopt.h \
	getopt1.c
liblalsupport_la_LIBADD = @LALSUPPORT_LIBADD@
liblalsupport_la_LDFLAGS = -version-info $(LIBVERSION_SUPPORT)

noinst_HEADERS = \
	PrintFrequencySeries_source.c \
	PrintTimeSeries_source.c \
	PrintVector_source.c \
	ReadFrequencySeries_source.c \
	ReadTimeSeries_source.c \
	StreamSequenceInput_source.c \
	StreamSequenceInputComplex_source.c \
	StreamVectorInput_source.c \
	StreamVectorSequenceInput_source.c

MAINTAINERCLEANFILES = \
	StreamGridInput.c \
	StreamGridOutput.c \
	StreamSeriesInput.c \
	StreamSeriesOutput.c

EXTRA_DIST = \
	LALReadGrid.m4 \
	LALReadSeries.m4 \
	LALWriteGrid.m4 \
	LALWriteSeries.m4 \
	StreamGridInput.m4 \
	StreamGridOutput.m4 \
	StreamSeriesInput.m4 \
	StreamSeriesOutput.m4

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

StreamSeriesInput.c: StreamSeriesInput.m4 LALReadSeries.m4
	cd $(srcdir) && @M4@ StreamSeriesInput.m4 > StreamSeriesInput.c

StreamSeriesOutput.c: StreamSeriesOutput.m4 LALWriteSeries.m4
	cd $(srcdir) && @M4@ StreamSeriesOutput.m4 > StreamSeriesOutput.c

StreamGridInput.c: StreamGridInput.m4 LALReadGrid.m4
	cd $(srcdir) && @M4@ StreamGridInput.m4 > StreamGridInput.c

StreamGridOutput.c: StreamGridOutput.m4 LALWriteGrid.m4
	cd $(srcdir) && @M4@ StreamGridOutput.m4 > StreamGridOutput.c

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

@DVIMKFILE@
