## Process this file with automake to produce Makefile.in
check_PROGRAMS = UnitsTest
UnitsTest_SOURCES = UnitsTest.c
UnitsTest_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = UnitsTest
MOSTLYCLEANFILES = PrintVector.*
