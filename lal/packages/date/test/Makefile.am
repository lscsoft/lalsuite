include $(top_srcdir)/gnuscripts/lalsuite_header_links.am

TEST_CPPFLAGS = -DTEST_DATA_DIR='"$(abs_srcdir)/"'

check_PROGRAMS = \
	TestJulianDay\
	TestLMST \
	TestLMST2 \
	TestDelay \
	TestGPStoUTC \
	TestUTCtoGPS \
	TestGPStoGMST1 \
	TestLeapSecs \
	TestGPStoFloat \
	TestIncrementGPS \
	TestStrToGPS

TESTS = \
	$(check_PROGRAMS)
