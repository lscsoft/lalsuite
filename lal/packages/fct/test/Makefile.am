## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALfctTest LALfctSample
LALfctTest_SOURCES = LALfctTest.c timingfunctions.c timingfunctions.h
LALfctTest_LDADD = $(top_builddir)/lib/liblal.la
LALfctSample_SOURCES = LALfctSample.c
LALfctSample_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = LALfctTest
MOSTLYCLEANFILES = output.dat
