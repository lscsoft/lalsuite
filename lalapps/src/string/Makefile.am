## Process this file with automake to produce Makefile.in
# $Id$
RINGSRC = ring.c
STRINGSRC = stringfind.c
if LALFRAME
if LALMETAIO
bin_PROGRAMS = lalapps_ring lalapps_stringfind
lalapps_ring_SOURCES = $(RINGSRC)
lalapps_stringfind_SOURCES = $(STRINGSRC)
TESTS = ring_test.sh
else
bin_PROGRAMS = lalapps_stringfind
lalapps_stringfind_SOURCES = $(STRINGSRC)
TESTS = 
endif
else
bin_PROGRAMS = lalapps_stringfind
lalapps_stringfind_SOURCES = $(STRINGSRC)
TESTS = 
endif
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
MOSTLYCLEANFILES = snr-* *.dat *.out *.xml
EXTRA_DIST = $(RINGSRC) $(STRINGSRC) \
             F-600000000-6-10.F F-600000060-6-10.F F-600000120-6-10.F \
             H-CAL_FAC-600000000-6000.gwf H-CAL_REF-600000000-64.gwf \
             H-600000000-600000180.cache H1-CAL-600000000-600006000.cache \
             response.asc filterpar.in ring_test.sh string_test.sh
