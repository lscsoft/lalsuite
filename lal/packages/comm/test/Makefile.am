## Process this file with automake to produce Makefile.in
if MPI
check_PROGRAMS  = SendRecvTest ExchangeTest
TESTS           = CommTest@MPI_TYPE@
else
check_PROGRAMS  =
TESTS           = NoTest
endif

EXTRA_DIST      = NoTest CommTestMPICH CommTestLAM SendRecvTest.c ExchangeTest.c

SendRecvTest_SOURCES = SendRecvTest.c
SendRecvTest_LDADD = $(top_builddir)/packages/comm/src/liblalmpi.la $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
SendRecvTest_CPPFLAGS = $(MPI_CPPFLAGS)
SendRecvTest_CFLAGS = $(MPI_CFLAGS)
SendRecvTest_LDFLAGS = $(MPI_LDFLAGS)

ExchangeTest_SOURCES = ExchangeTest.c
ExchangeTest_LDADD = $(top_builddir)/packages/comm/src/liblalmpi.la $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
ExchangeTest_CPPFLAGS = $(MPI_CPPFLAGS)
ExchangeTest_CFLAGS = $(MPI_CFLAGS)
ExchangeTest_LDFLAGS = $(MPI_LDFLAGS)


INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep PI* Slave*.out Slave*.err machines killscript lamhosts CommTestSchema
CLEANFILES = SendRecvTest ExchangeTest
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
