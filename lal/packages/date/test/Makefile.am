AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include -I$(srcdir)/../src
LDADD = \
	$(top_builddir)/lib/lal/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

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

TESTS_ENVIRONMENT = 1>/dev/null
