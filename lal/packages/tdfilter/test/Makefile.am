## Process this file with automake to produce Makefile.in
check_PROGRAMS = BandPassTest IIRFilterTest
BandPassTest_SOURCES = BandPassTest.c
BandPassTest_LDADD = $(top_srcdir)/lib/liblal.la
IIRFilterTest_SOURCES = IIRFilterTest.c
IIRFilterTest_LDADD = $(top_srcdir)/lib/liblal.la
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
TESTS = BandPassTest IIRFilterTest
MOSTLYCLEANFILES = out.dat
