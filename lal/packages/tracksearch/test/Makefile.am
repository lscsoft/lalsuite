## Process this file with automake to produce Makefile.in
check_PROGRAMS = TrackSearchTest
TrackSearchTest_SOURCES = TrackSearchTest.c
TrackSearchTest_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = TrackSearchTestScript.sh
EXTRA_DIST = a.pgm TrackSearchTestScript.sh
