## Process this file with automake to produce Makefile.in
noinst_PROGRAMS = BandPassTest IIRFilterTest
BandPassTest_SOURCES = BandPassTest.c
BandPassTest_LDADD = $(top_srcdir)/lib/liblal.a
IIRFilterTest_SOURCES = IIRFilterTest.c
IIRFilterTest_LDADD = $(top_srcdir)/lib/liblal.a
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
TESTS = $(PROGRAMS)
MOSTLYCLEANFILES = out.dat
