## Process this file with automake to produce Makefile.in
COMMSOURCES = SendRecv.c Exchange.c

if MPI
lib_LTLIBRARIES = liblalmpi.la
liblalmpi_la_SOURCES = $(COMMSOURCES)
liblalmpi_la_CPPFLAGS = $(MPI_CPPFLAGS)
liblalmpi_la_CFLAGS = $(MPI_CFLAGS)
liblalmpi_la_LIBADD = $(MPI_LDFLAGS)
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
	cd $(srcdir) && @M4@ SendRecv.m4 > SendRecv.c

Exchange.c: Exchange.m4 ExchangeObject.m4
	cd $(srcdir) && @M4@ Exchange.m4 > Exchange.c
@DVIMKFILE@
