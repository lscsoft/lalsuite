## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libsupport.la
libsupport_la_SOURCES = LALVersion.c LALMalloc.c LALError.c PrintVector.c \
                        PrintFrequencySeries.c PrintTimeSeries.c
MAINTAINERCLEANFILES = PrintVector.c \
                        PrintFrequencySeries.c PrintTimeSeries.c
EXTRA_DIST = PrintVector.m4 LALPrintVector.m4 \
             PrintTimeSeries.m4 LALPrintTimeSeries.m4 \
             PrintFrequencySeries.m4 LALPrintFrequencySeries.m4
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
PrintVector.c: PrintVector.m4 LALPrintVector.m4
	m4 PrintVector.m4 > PrintVector.c
PrintTimeSeries.c: PrintTimeSeries.m4 LALPrintTimeSeries.m4
	m4 PrintTimeSeries.m4 > PrintTimeSeries.c
PrintFrequencySeries.c: PrintFrequencySeries.m4 LALPrintFrequencySeries.m4
	m4 PrintFrequencySeries.m4 > PrintFrequencySeries.c
