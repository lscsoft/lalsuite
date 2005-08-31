## Process this file with automake to produce Makefile.in
#
# $Id$
#
check_PROGRAMS = TestJulianDay TestLMST TestLMST2 TestDelay TestDateString \
                 TestGPStoUTC TestUTCtoGPS TestGPStoGMST1 TestLeapSecs \
		 TestGPStoFloat TestIncrementGPS TestStrToGPS
TestJulianDay_SOURCES = TestJulianDay.c
TestJulianDay_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestLMST_SOURCES = TestLMST.c
TestLMST_LDADD = $(top_builddir)/lib/liblal.la
TestLMST2_SOURCES = TestLMST2.c
TestLMST2_LDADD = $(top_builddir)/lib/liblal.la
TestDelay_SOURCES = TestDelay.c
TestDelay_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestDateString_SOURCES = TestDateString.c
TestDateString_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestGPStoUTC_SOURCES = TestGPStoUTC.c
TestGPStoUTC_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestUTCtoGPS_SOURCES = TestUTCtoGPS.c
TestUTCtoGPS_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestGPStoGMST1_SOURCES = TestGPStoGMST1.c
TestGPStoGMST1_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestLeapSecs_SOURCES = TestLeapSecs.c
TestLeapSecs_LDADD   = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestGPStoFloat_SOURCES = TestGPStoFloat.c
TestGPStoFloat_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
#TestGPStoMST_SOURCES = TestGPStoMST.c
#TestGPStoMST_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestIncrementGPS_SOURCES = TestIncrementGPS.c
TestIncrementGPS_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestStrToGPS_SOURCES = TestStrToGPS.c
TestStrToGPS_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = TestJulianDay TestLMST TestDelay TestDateString TestGPStoUTC \
        TestUTCtoGPS TestLeapSecs TestGPStoFloat TestIncrementGPS TestStrToGPS
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
