## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALInspiralTest
LALInspiralTest_SOURCES = LALInspiralTest.c
LALInspiralTest_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = LALInspiralTest
CLEANFILES = *.dat
