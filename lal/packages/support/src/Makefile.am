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
	VeryBasicXMLOutput.c \
	getopt.c \
	getopt.h \
	getopt1.c
liblalsupport_la_LIBADD = @LALSUPPORT_LIBADD@
liblalsupport_la_LDFLAGS = -version-info $(LIBVERSION_SUPPORT)

MAINTAINERCLEANFILES = \
	PrintFrequencySeries.c \
	PrintTimeSeries.c \
	PrintVector.c \
	ReadFrequencySeries.c \
	ReadTimeSeries.c \
	StreamGridInput.c \
	StreamGridOutput.c \
	StreamSequenceInput.c \
	StreamSeriesInput.c \
	StreamSeriesOutput.c \
	StreamVectorInput.c \
	StreamVectorSequenceInput.c

EXTRA_DIST = \
	LALPrintFrequencySeries.m4 \
	LALPrintTimeSeries.m4 \
	LALPrintVector.m4 \
	LALReadComplexSequence.m4 \
	LALReadFrequencySeries.m4 \
	LALReadGrid.m4 \
	LALReadSequence.m4 \
	LALReadSeries.m4 \
	LALReadTimeSeries.m4 \
	LALReadVector.m4 \
	LALReadVectorSequence.m4 \
	LALWriteGrid.m4 \
	LALWriteSeries.m4 \
	PrintFrequencySeries.m4 \
	PrintTimeSeries.m4 \
	PrintVector.m4 \
	ReadFrequencySeries.m4 \
	ReadTimeSeries.m4 \
	StreamGridInput.m4 \
	StreamGridOutput.m4 \
	StreamSeriesInput.m4 \
	StreamSeriesOutput.m4 \
	StreamSequenceInput.m4 \
	StreamVectorInput.m4 \
	StreamVectorSequenceInput.m4 \
	snprintf_LICENSE.txt \
	snprintf_README

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

PrintVector.c: PrintVector.m4 LALPrintVector.m4
	cd $(srcdir) && @M4@ PrintVector.m4 > PrintVector.c

PrintTimeSeries.c: PrintTimeSeries.m4 LALPrintTimeSeries.m4
	cd $(srcdir) && @M4@ PrintTimeSeries.m4 > PrintTimeSeries.c

PrintFrequencySeries.c: PrintFrequencySeries.m4 LALPrintFrequencySeries.m4
	cd $(srcdir) && @M4@ PrintFrequencySeries.m4 > PrintFrequencySeries.c

ReadFrequencySeries.c: ReadFrequencySeries.m4 LALReadFrequencySeries.m4
	cd $(srcdir) && @M4@ ReadFrequencySeries.m4 > ReadFrequencySeries.c

ReadTimeSeries.c: ReadTimeSeries.m4 LALReadTimeSeries.m4
	cd $(srcdir) && @M4@ ReadTimeSeries.m4 > ReadTimeSeries.c

StreamVectorInput.c: StreamVectorInput.m4 LALReadVector.m4
	cd $(srcdir) && @M4@ StreamVectorInput.m4 > StreamVectorInput.c

StreamVectorSequenceInput.c: StreamVectorSequenceInput.m4 LALReadVectorSequence.m4
	cd $(srcdir) && @M4@ StreamVectorSequenceInput.m4 > StreamVectorSequenceInput.c

StreamSequenceInput.c: StreamSequenceInput.m4 LALReadSequence.m4 LALReadComplexSequence.m4
	cd $(srcdir) && @M4@ StreamSequenceInput.m4 > StreamSequenceInput.c

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
