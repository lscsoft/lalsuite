## Process this file with automake to produce Makefile.in
check_PROGRAMS = TrackSearchTest
TrackSearchTest_SOURCES = TrackSearchTest.c
TrackSearchTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = TrackSearchTest
TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null LAL_DATA_PATH=$(srcdir)
noinst_DATA = a.pgm
EXTRA_DIST = $(DATA)
@DVIMKFILE@
