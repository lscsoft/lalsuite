## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libcomm.la
COMM_SOURCES = SendRecv.c Exchange.c

if MPI
libcomm_la_SOURCES = $(COMM_SOURCES)
EXTRA_TEXSOURCES   =
else
libcomm_la_SOURCES = 
EXTRA_TEXSOURCES   = $(COMM_SOURCES)
endif

EXTRA_DIST = $(COMM_SOURCES) \
  SendRecv.m4 SendRecvAtomic.m4 SendRecvVector.m4 SendRecvSeries.m4 \
  Exchange.m4 ExchangeObject.m4
MAINTAINERCLEANFILES = $(COMM_SOURCES)

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

SendRecv.c: SendRecv.m4 SendRecvAtomic.m4 SendRecvVector.m4 SendRecvSeries.m4
	m4 SendRecv.m4 > SendRecv.c

Exchange.c: Exchange.m4 ExchangeObject.m4
	m4 Exchange.m4 > Exchange.c
