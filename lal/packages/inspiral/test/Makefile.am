## Process this file with automake to produce Makefile.in
check_PROGRAMS = InspiralTest
InspiralTest_SOURCES = InspiralTest.c
InspiralTest_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = InspiralTest
CLEANFILES = *.dat
