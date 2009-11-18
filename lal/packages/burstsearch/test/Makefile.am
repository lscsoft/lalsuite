## Process this file with automake to produce Makefile.in
check_PROGRAMS = ThresholdsTest

ThresholdsTest_SOURCES = ThresholdsTest.c
ThresholdsTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = ThresholdsTest
TESTS_ENVIRONMENT = 1>/dev/null
MOSTLYCLEANFILES = PrintVector.*
