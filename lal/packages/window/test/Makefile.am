## Process this file with automake to produce Makefile.in
check_PROGRAMS = WindowTest
WindowTest_SOURCES = WindowTest.c
WindowTest_LDADD = $(top_builddir)/lib/liblal.la -lm
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = WindowTest
TESTS_ENVIRONMENT = 1>/dev/null
