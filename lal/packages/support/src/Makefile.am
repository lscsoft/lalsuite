## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libsupport.la
lib_LTLIBRARIES = liblalsupport.la

libsupport_la_SOURCES = LALStd.c LALVersion.c LALMalloc.c LALError.c \
  StringToken.c StringConvert.c

liblalsupport_la_SOURCES = FileIO.c PrintVector.c PrintFrequencySeries.c \
  PrintTimeSeries.c ReadFrequencySeries.c ReadTimeSeries.c \
  StreamVectorInput.c StreamVectorSequenceInput.c StreamSequenceInput.c \
  StreamSeriesInput.c StreamSeriesOutput.c \
  LALInitBarycenter.c LALXMGRInterface.c

MAINTAINERCLEANFILES = PrintVector.c PrintFrequencySeries.c PrintTimeSeries.c \
  ReadFrequencySeries.c ReadTimeSeries.c \
  StreamVectorInput.c StreamVectorSequenceInput.c StreamSequenceInput.c \
  StreamSeriesInput.c StreamSeriesOutput.c

EXTRA_DIST = PrintVector.m4 LALPrintVector.m4 \
             PrintTimeSeries.m4 LALPrintTimeSeries.m4 \
	     PrintFrequencySeries.m4 LALPrintFrequencySeries.m4\
	     ReadFrequencySeries.m4 LALReadFrequencySeries.m4\
             snprintf_LICENSE.txt snprintf_README \
             ReadTimeSeries.m4 LALReadTimeSeries.m4 \
             StreamVectorInput.m4 LALReadVector.m4 \
             StreamVectorSequenceInput.m4 LALReadVectorSequence.m4 \
             StreamSequenceInput.m4 LALReadSequence.m4 LALReadComplexSequence.m4 \
             StreamSeriesInput.m4 LALReadSeries.m4 \
             StreamSeriesOutput.m4 LALWriteSeries.m4

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

PrintVector.c: PrintVector.m4 LALPrintVector.m4
	cd $(srcdir) && m4 PrintVector.m4 > PrintVector.c

PrintTimeSeries.c: PrintTimeSeries.m4 LALPrintTimeSeries.m4
	cd $(srcdir) && m4 PrintTimeSeries.m4 > PrintTimeSeries.c

PrintFrequencySeries.c: PrintFrequencySeries.m4 LALPrintFrequencySeries.m4
	cd $(srcdir) && m4 PrintFrequencySeries.m4 > PrintFrequencySeries.c

ReadFrequencySeries.c: ReadFrequencySeries.m4 LALReadFrequencySeries.m4
	cd $(srcdir) && m4 ReadFrequencySeries.m4 > ReadFrequencySeries.c

ReadTimeSeries.c: ReadTimeSeries.m4 LALReadTimeSeries.m4
	cd $(srcdir) && m4 ReadTimeSeries.m4 > ReadTimeSeries.c

StreamVectorInput.c: StreamVectorInput.m4 LALReadVector.m4
	cd $(srcdir) && m4 StreamVectorInput.m4 > StreamVectorInput.c

StreamVectorSequenceInput.c: StreamVectorSequenceInput.m4 LALReadVectorSequence.m4
	cd $(srcdir) && m4 StreamVectorSequenceInput.m4 > StreamVectorSequenceInput.c

StreamSequenceInput.c: StreamSequenceInput.m4 LALReadSequence.m4 LALReadComplexSequence.m4
	cd $(srcdir) && m4 StreamSequenceInput.m4 > StreamSequenceInput.c

StreamSeriesInput.c: StreamSeriesInput.m4 LALReadSeries.m4
	cd $(srcdir) && m4 StreamSeriesInput.m4 > StreamSeriesInput.c

StreamSeriesOutput.c: StreamSeriesOutput.m4 LALWriteSeries.m4
	cd $(srcdir) && m4 StreamSeriesOutput.m4 > StreamSeriesOutput.c
