## Process this file with automake to produce Makefile.in
# $Id$
bin_PROGRAMS = lalapps_ring
lalapps_ring_SOURCES = ring.c
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
MOSTLYCLEANFILES = snr-* *.dat *.out
EXTRA_DIST = F-600000000-6-10.F F-600000060-6-10.F F-600000120-6-10.F \
  response.asc filterpar.in ring_test.sh
TESTS = ring_test.sh
