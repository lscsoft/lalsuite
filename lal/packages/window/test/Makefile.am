## Process this file with automake to produce Makefile.in
check_PROGRAMS = WindowTest
WindowTest_SOURCES = WindowTest.c
WindowTest_LDADD = $(top_srcdir)/lib/liblal.la -lm
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
TESTS = WindowTest
MOSTLYCLEANFILES = PrintVector.*
