## Process this file with automake to produce Makefile.in
check_PROGRAMS = ComplexFFTTest RealFFTTest TimeFreqFFTTest AverageSpectrumTest AvgSpecTest
TimeFreqFFTTest_SOURCES = TimeFreqFFTTest.c
ComplexFFTTest_SOURCES = ComplexFFTTest.c
RealFFTTest_SOURCES = RealFFTTest.c
AverageSpectrumTest_SOURCES = AverageSpectrumTest.c
AvgSpecTest_SOURCES = AvgSpecTest.c
RealFFTTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
ComplexFFTTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TimeFreqFFTTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
AverageSpectrumTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
AvgSpecTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = ComplexFFTTest RealFFTTest TimeFreqFFTTest AverageSpectrumTest AvgSpecTest
TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
MOSTLYCLEANFILES = .dvi-dep *.out out*.dat
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
