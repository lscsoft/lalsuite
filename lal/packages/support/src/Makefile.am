## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libsupport.la
lib_LTLIBRARIES = liblalsupport.la

libsupport_la_SOURCES = LALStd.c LALVersion.c LALMalloc.c LALError.c

liblalsupport_la_SOURCES = FileIO.c PrintVector.c PrintFrequencySeries.c \
  PrintTimeSeries.c ReadFrequencySeries.c ReadTimeSeries.c \
  StreamVectorInput.c StreamVectorSequenceInput.c LALInitBarycenter.c \
  LALXMGRInterface.c $(top_srcdir)/packages/findchirp/src/FindChirpView.c 

MAINTAINERCLEANFILES = PrintVector.c PrintFrequencySeries.c PrintTimeSeries.c \
  ReadFrequencySeries.c ReadTimeSeries.c

EXTRA_DIST = PrintVector.m4 LALPrintVector.m4 \
             PrintTimeSeries.m4 LALPrintTimeSeries.m4 \
	     PrintFrequencySeries.m4 LALPrintFrequencySeries.m4\
	     ReadFrequencySeries.m4 LALReadFrequencySeries.m4\
             snprintf_LICENSE.txt snprintf_README \
             ReadTimeSeries.m4 LALReadTimeSeries.m4

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

PrintVector.c: PrintVector.m4 LALPrintVector.m4
	m4 PrintVector.m4 > PrintVector.c

PrintTimeSeries.c: PrintTimeSeries.m4 LALPrintTimeSeries.m4
	m4 PrintTimeSeries.m4 > PrintTimeSeries.c

PrintFrequencySeries.c: PrintFrequencySeries.m4 LALPrintFrequencySeries.m4
	m4 PrintFrequencySeries.m4 > PrintFrequencySeries.c

ReadFrequencySeries.c: ReadFrequencySeries.m4 LALReadFrequencySeries.m4
	m4 ReadFrequencySeries.m4 > ReadFrequencySeries.c

ReadTimeSeries.c: ReadTimeSeries.m4 LALReadTimeSeries.m4
	m4 ReadTimeSeries.m4 > ReadTimeSeries.c

