## Process this file with automake to produce Makefile.in
noinst_PROGRAMS = RealFFTTest RealPowerSpectrumTest
RealFFTTest_SOURCES = RealFFTTest.c
RealPowerSpectrumTest_SOURCES = RealPowerSpectrumTest.c
RealFFTTest_LDADD = $(top_srcdir)/lib/liblal.a
RealPowerSpectrumTest_LDADD = $(top_srcdir)/lib/liblal.a
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
TESTS = $(PROGRAMS)
MOSTLYCLEANFILES = RealPowerSpectrumTest.out
