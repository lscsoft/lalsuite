AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include -I$(srcdir)/../src @AM_CPPFLAGS@
LDADD = \
	$(top_builddir)/lib/lal/liblal.la \
	$(top_builddir)/lib/lalsupport/src/liblalsupport.la

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
