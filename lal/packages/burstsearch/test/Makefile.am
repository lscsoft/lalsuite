## Process this file with automake to produce Makefile.in
check_PROGRAMS = ThresholdsTest TFTransformTest ExcessPowerTest

ThresholdsTest_SOURCES = ThresholdsTest.c
ThresholdsTest_LDADD = $(top_srcdir)/lib/liblal.la

TFTransformTest_SOURCES = TFTransformTest.c
TFTransformTest_LDADD = $(top_srcdir)/lib/liblal.la

ExcessPowerTest_SOURCES = ExcessPowerTest.c
ExcessPowerTest_LDADD = $(top_srcdir)/lib/liblal.la


INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
TESTS = ThresholdsTest TFTransformTest ExcessPowerTest
MOSTLYCLEANFILES = PrintVector.*
