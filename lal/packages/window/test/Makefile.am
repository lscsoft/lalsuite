## Process this file with automake to produce Makefile.in
noinst_PROGRAMS = WindowTest
WindowTest_SOURCES = WindowTest.c
WindowTest_LDADD = $(top_srcdir)/lib/liblal.a -lm
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
TESTS = $(PROGRAMS)
MOSTLYCLEANFILES = PrintVector.*
