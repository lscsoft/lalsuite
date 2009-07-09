## Process this file with automake to produce Makefile.in
check_PROGRAMS = ThresholdsTest SkymapTest

ThresholdsTest_SOURCES = ThresholdsTest.c
ThresholdsTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

SkymapTest_SOURCES = SkymapTest.c
SkymapTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = ThresholdsTest SkymapTest
TESTS_ENVIRONMENT = 1>/dev/null
MOSTLYCLEANFILES = PrintVector.*
