## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALfctTest LALfctSample LALFCTInterfaceTest
LALfctTest_SOURCES = LALfctTest.c timingfunctions.c timingfunctions.h
LALfctTest_LDADD = $(top_builddir)/lib/liblal.la
LALfctSample_SOURCES = LALfctSample.c
LALfctSample_LDADD = $(top_builddir)/lib/liblal.la
LALFCTInterfaceTest_SOURCES = LALFCTInterfaceTest.c
LALFCTInterfaceTest_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
#TESTS = LALfctTest LALFCTInterfaceTest
TESTS = LALfctTest
MOSTLYCLEANFILES = output.dat
