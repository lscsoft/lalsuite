## Process this file with automake to produce Makefile.in
#check_PROGRAMS = ThresholdsTest TFTransformTest ExcessPowerTest
check_PROGRAMS = ThresholdsTest TFTransformTest

ThresholdsTest_SOURCES = ThresholdsTest.c
ThresholdsTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la

TFTransformTest_SOURCES = TFTransformTest.c
TFTransformTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la

#ExcessPowerTest_SOURCES = ExcessPowerTest.c
#ExcessPowerTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la


INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
#TESTS = ThresholdsTest TFTransformTest ExcessPowerTest
TESTS = ThresholdsTest TFTransformTest
MOSTLYCLEANFILES = .dvi-dep PrintVector.*
DISTCLEANFILES = *.tex
EXTRA_DIST = ExcessPowerTest.c
EXTRA_TEXSOURCES = ExcessPowerTest.c
dvi-local: .dvi-dep
