## Process this file with automake to produce Makefile.in
check_PROGRAMS = ComplexFFTTest RealFFTTest TimeFreqFFTTest
TimeFreqFFTTest_SOURCES = TimeFreqFFTTest.c
ComplexFFTTest_SOURCES = ComplexFFTTest.c
RealFFTTest_SOURCES = RealFFTTest.c
RealFFTTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
ComplexFFTTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TimeFreqFFTTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = ComplexFFTTest RealFFTTest TimeFreqFFTTest
MOSTLYCLEANFILES = .dvi-dep *.out
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
