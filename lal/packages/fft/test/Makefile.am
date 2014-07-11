include $(top_srcdir)/gnuscripts/lallibs_header_links.am

TEST_CPPFLAGS = -DTEST_DATA_DIR='"$(abs_srcdir)/"'

check_PROGRAMS = \
	ComplexFFTTest \
	RealFFTTest \
	TimeFreqFFTTest \
	AverageSpectrumTest \
	AvgSpecTest

TESTS = \
	$(check_PROGRAMS)

MOSTLYCLEANFILES = *.out out*.dat
