## Process this file with automake to produce Makefile.in
noinst_PROGRAMS = ThresholdsTest TFTransformTest ExcessPowerTest

ThresholdsTest_SOURCES = ThresholdsTest.c
ThresholdsTest_LDADD = $(top_srcdir)/lib/liblal.a

TFTransformTest_SOURCES = TFTransformTest.c
TFTransformTest_LDADD = $(top_srcdir)/lib/liblal.a

ExcessPowerTest_SOURCES = ExcessPowerTest.c
ExcessPowerTest_LDADD = $(top_srcdir)/lib/liblal.a


INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
TESTS = $(PROGRAMS)
MOSTLYCLEANFILES = PrintVector.*
