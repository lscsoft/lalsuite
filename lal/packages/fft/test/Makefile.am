## Process this file with automake to produce Makefile.in
check_PROGRAMS = ComplexFFTTest RealFFTTest RealPowerSpectrumTest
ComplexFFTTest_SOURCES = ComplexFFTTest.c
RealFFTTest_SOURCES = RealFFTTest.c
RealPowerSpectrumTest_SOURCES = RealPowerSpectrumTest.c
RealFFTTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
ComplexFFTTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
RealPowerSpectrumTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = ComplexFFTTest RealFFTTest RealPowerSpectrumTest
MOSTLYCLEANFILES = .dvi-dep RealPowerSpectrumTest.out
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
