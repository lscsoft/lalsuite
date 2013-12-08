AM_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_builddir)/include -I$(srcdir)/../include
LDADD = \
	$(top_builddir)/lib/lal/liblal.la \
	$(top_builddir)/lib/lalsupport/src/liblalsupport.la

check_PROGRAMS = \
	ComplexFFTTest \
	RealFFTTest \
	TimeFreqFFTTest \
	AverageSpectrumTest \
	AvgSpecTest

TESTS = \
	$(check_PROGRAMS)

MOSTLYCLEANFILES = *.out out*.dat
