## Process this file with automake to produce Makefile.in
check_PROGRAMS = TrackSearchTest
TrackSearchTest_SOURCES = TrackSearchTest.c
TrackSearchTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = TrackSearchTest
TESTS_ENVIRONMENT = LAL_DATA_PATH=$(srcdir)
noinst_DATA = a.pgm
EXTRA_DIST = $(DATA)
