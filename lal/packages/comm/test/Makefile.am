## Process this file with automake to produce Makefile.in
if MPI
check_PROGRAMS  = SendRecvTest ExchangeTest
TESTS           = CommTest@MPITYPE@
else
check_PROGRAMS  =
TESTS           = NoTest
endif

EXTRA_DIST      = NoTest CommTestMPICH CommTestLAM SendRecvTest.c ExchangeTest.c

SendRecvTest_SOURCES = SendRecvTest.c
SendRecvTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
ExchangeTest_SOURCES = ExchangeTest.c
ExchangeTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep PI* Slave*.out Slave*.err machines killscript lamhosts CommTestSchema
CLEANFILES = SendRecvTest ExchangeTest
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
