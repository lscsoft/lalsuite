## Process this file with automake to produce Makefile.in
#check_PROGRAMS = TestJulianDay TestUtoGPS TestLMST TestDelay
check_PROGRAMS = TestJulianDay TestUtoGPS TestDelay
TestJulianDay_SOURCES = TestJulianDay.c
TestJulianDay_LDADD = $(top_builddir)/lib/liblal.la
TestUtoGPS_SOURCES = TestUtoGPS.c
TestUtoGPS_LDADD = $(top_builddir)/lib/liblal.la
#TestLMST_SOURCES = TestLMST.c
#TestLMST_LDADD = $(top_builddir)/lib/liblal.la
TestDelay_SOURCES = TestDelay.c
TestDelay_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
#TESTS = TestJulianDay TestUtoGPS TestLMST TestDelay
TESTS = TestJulianDay TestUtoGPS TestDelay
