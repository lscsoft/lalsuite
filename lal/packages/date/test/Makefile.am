## Process this file with automake to produce Makefile.in
check_PROGRAMS = TestJulianDay TestUtoGPS TestLMST TestLMST2 TestDelay TestDateString
TestJulianDay_SOURCES = TestJulianDay.c
TestJulianDay_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
TestUtoGPS_SOURCES = TestUtoGPS.c
TestUtoGPS_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
TestLMST_SOURCES = TestLMST.c
TestLMST_LDADD = $(top_builddir)/lib/liblal.la
TestLMST2_SOURCES = TestLMST2.c
TestLMST2_LDADD = $(top_builddir)/lib/liblal.la
TestDelay_SOURCES = TestDelay.c
TestDelay_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
TestDateString_SOURCES = TestDateString.c
TestDateString_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
#TestGPStoMST_SOURCES = TestGPStoMST.c
#TestGPStoMST_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
#TESTS = TestJulianDay TestUtoGPS TestLMST TestDelay TestGPStoMST
TESTS = TestJulianDay TestUtoGPS TestLMST TestDelay TestDateString
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
