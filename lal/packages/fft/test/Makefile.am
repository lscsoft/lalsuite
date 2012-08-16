INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

check_PROGRAMS = \
	ComplexFFTTest \
	RealFFTTest \
	TimeFreqFFTTest \
	AverageSpectrumTest \
	AvgSpecTest

TESTS = \
	$(check_PROGRAMS)

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null

MOSTLYCLEANFILES = *.out out*.dat
