## Process this file with automake to produce Makefile.in
check_PROGRAMS = TFClustersTest1 TFClustersTest2 
TFClustersTest1_SOURCES = TFClustersTest1.c
TFClustersTest1_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la

TFClustersTest2_SOURCES = TFClustersTest2.c
TFClustersTest2_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = TFClustersTest1 TFClustersTest2 
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
