## Process this file with automake to produce Makefile.in
noinst_PROGRAMS = TestJulianDay TestUtoGPS TestLMST
TestJulianDay_SOURCES = TestJulianDay.c
TestJulianDay_LDADD = $(top_srcdir)/lib/liblal.a
TestUtoGPS_SOURCES = TestUtoGPS.c
TestUtoGPS_LDADD = $(top_srcdir)/lib/liblal.a
TestLMST_SOURCES = TestLMST.c
TestLMST_LDADD = $(top_srcdir)/lib/liblal.a
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
TESTS = $(PROGRAMS)
