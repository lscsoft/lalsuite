## Process this file with automake to produce Makefile.in
#check_PROGRAMS = TFCThresholdsTest1 TFClustersTest1 TFClustersTest2 TFClustersTest3 TFClustersTest4 TFCWaveletTest1
check_PROGRAMS = TFCThresholdsTest1 TFClustersTest1 TFClustersTest2 TFClustersTest3 TFClustersTest4

TFClustersTest1_SOURCES = TFClustersTest1.c
TFClustersTest1_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

TFClustersTest2_SOURCES = TFClustersTest2.c
TFClustersTest2_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

TFClustersTest3_SOURCES = TFClustersTest3.c
TFClustersTest3_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

TFClustersTest4_SOURCES = TFClustersTest4.c
TFClustersTest4_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

TFCThresholdsTest1_SOURCES = TFCThresholdsTest1.c
TFCThresholdsTest1_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

#TFCWaveletTest1_SOURCES = TFCWaveletTest1.c
#TFCWaveletTest1_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

EXTRA_DIST = TFCWaveletTest1.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
#TESTS = TFCThresholdsTest1 TFClustersTest1 TFClustersTest2 TFClustersTest3 TFClustersTest4 TFCWaveletTest1
TESTS = TFCThresholdsTest1 TFClustersTest1 TFClustersTest2 TFClustersTest3 TFClustersTest4
TESTS_ENVIRONMENT = 1>/dev/null
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
