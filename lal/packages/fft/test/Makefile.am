## Process this file with automake to produce Makefile.in
check_PROGRAMS = ComplexFFTTest RealFFTTest RealPowerSpectrumTest
ComplexFFTTest_SOURCES = ComplexFFTTest.c
RealFFTTest_SOURCES = RealFFTTest.c
RealPowerSpectrumTest_SOURCES = RealPowerSpectrumTest.c
RealFFTTest_LDADD = $(top_srcdir)/lib/liblal.la
ComplexFFTTest_LDADD = $(top_srcdir)/lib/liblal.la
RealPowerSpectrumTest_LDADD = $(top_srcdir)/lib/liblal.la
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
TESTS = ComplexFFTTest RealFFTTest RealPowerSpectrumTest
MOSTLYCLEANFILES = RealPowerSpectrumTest.out
