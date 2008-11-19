## Process this file with automake to produce Makefile.in
## $Id$
LALMETAIOSOURCES = CreateMetaTableDir.c LIGOLwXMLRead.c LIGOLwXMLRingdownRead.c
LALXMLSOURCES = LALXML.c

if METAIO
LALMETAIOLIBRARIES = liblalmetaio.la
else
LALMETAIOLIBRARIES = 
endif

if XML
LALXMLLIBRARIES = liblalxml.la
else
LALXMLLIBRARIES = 
endif

lib_LTLIBRARIES = $(LALMETAIOLIBRARIES) $(LALXMLLIBRARIES) liblalsupport.la

if METAIO
liblalmetaio_la_SOURCES = $(LALMETAIOSOURCES)
else
EXTRA_TEXSOURCES = $(LALMETAIOSOURCES)
endif

if XML
liblalxml_la_SOURCES = $(LALXMLSOURCES)
endif

liblalsupport_la_SOURCES = FileIO.c PrintVector.c PrintFrequencySeries.c \
  PrintTimeSeries.c ReadFrequencySeries.c ReadTimeSeries.c \
  ReadNoiseSpectrum.c \
  StreamVectorInput.c StreamVectorSequenceInput.c StreamSequenceInput.c \
  StreamSeriesInput.c StreamSeriesOutput.c \
  StreamGridInput.c StreamGridOutput.c \
  LALCache.c \
  LALInitBarycenter.c LALXMGRInterface.c \
  LIGOLwXML.c LIGOLwXMLArray.c SFTfileIO.c \
  LogPrintf.c UserInput.c ConfigFile.c getopt.c getopt.h getopt1.c \
  LALMath3DPlot.c LALMathNDPlot.c SegmentsIO.c \
  Audio.c NRWaveIO.c	\
  DopplerScan.c DopplerFullScan.c DopplerLatticeCovering.c \
  VeryBasicXMLOutput.c ReadFiltersFile.c

liblalsupport_la_LIBADD = @LALSUPPORT_LIBADD@


MAINTAINERCLEANFILES = PrintVector.c PrintFrequencySeries.c PrintTimeSeries.c \
  ReadFrequencySeries.c ReadTimeSeries.c \
  StreamVectorInput.c StreamVectorSequenceInput.c StreamSequenceInput.c \
  StreamSeriesInput.c StreamSeriesOutput.c \
  StreamGridInput.c StreamGridOutput.c

EXTRA_DIST = $(LALMETAIOSOURCES) $(LALXMLSOURCES)\
             PrintVector.m4 LALPrintVector.m4 \
             PrintTimeSeries.m4 LALPrintTimeSeries.m4 \
	     PrintFrequencySeries.m4 LALPrintFrequencySeries.m4\
	     ReadFrequencySeries.m4 LALReadFrequencySeries.m4\
             snprintf_LICENSE.txt snprintf_README \
             ReadTimeSeries.m4 LALReadTimeSeries.m4 \
             StreamVectorInput.m4 LALReadVector.m4 \
             StreamVectorSequenceInput.m4 LALReadVectorSequence.m4 \
             StreamSequenceInput.m4 LALReadSequence.m4 LALReadComplexSequence.m4 \
             StreamSeriesInput.m4 LALReadSeries.m4 \
             StreamSeriesOutput.m4 LALWriteSeries.m4 \
             StreamGridInput.m4 LALReadGrid.m4 \
             StreamGridOutput.m4 LALWriteGrid.m4

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

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
@DVIMKFILE@
