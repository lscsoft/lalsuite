## Process this file with automake to produce Makefile.in
COMMSOURCES = SendRecv.c Exchange.c

if MPI
noinst_LTLIBRARIES = libcomm.la
libcomm_la_SOURCES = $(COMMSOURCES)
EXTRA_TEXSOURCES   =
else
EXTRA_TEXSOURCES   = $(COMMSOURCES)
endif

EXTRA_DIST = $(COMMSOURCES) \
  SendRecv.m4 SendRecvAtomic.m4 SendRecvVector.m4 SendRecvSeries.m4 \
  Exchange.m4 ExchangeObject.m4
MAINTAINERCLEANFILES = $(COMMSOURCES)

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

SendRecv.c: SendRecv.m4 SendRecvAtomic.m4 SendRecvVector.m4 SendRecvSeries.m4
	cd $(srcdir) && m4 SendRecv.m4 > SendRecv.c

Exchange.c: Exchange.m4 ExchangeObject.m4
	cd $(srcdir) && m4 Exchange.m4 > Exchange.c
