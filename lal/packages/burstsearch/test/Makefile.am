## Process this file with automake to produce Makefile.in
check_PROGRAMS = SlopeETGTest TFClustersETGTest StdBurstOutput ThresholdsTest TFTransformTest

SlopeETGTest_SOURCES = SlopeETGTest.c
SlopeETGTest_LDADD = $(top_builddir)/lib/liblal.la 

TFClustersETGTest_SOURCES = TFClustersETGTest.c
TFClustersETGTest_LDADD = $(top_builddir)/lib/liblal.la 

StdBurstOutput_SOURCES = StdBurstOutput.c
StdBurstOutput_LDADD = $(top_builddir)/lib/liblal.la

ThresholdsTest_SOURCES = ThresholdsTest.c
ThresholdsTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

TFTransformTest_SOURCES = TFTransformTest.c
TFTransformTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

ExcessPowerTest_SOURCES = ExcessPowerTest.c
ExcessPowerTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la


INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = SlopeETGTest StdBurstOutput TFClustersETGTest ThresholdsTest TFTransformTest
MOSTLYCLEANFILES = .dvi-dep PrintVector.*
DISTCLEANFILES = *.tex
EXTRA_DIST = ExcessPowerTest.c
EXTRA_TEXSOURCES = ExcessPowerTest.c
dvi-local: .dvi-dep
@DVIMKFILE@
