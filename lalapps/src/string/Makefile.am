## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
if LALFRAME
if LALMETAIO
RINGPROG = lalapps_ring
else
RINGPROG =
endif
else
RINGPROG =
endif
else
RINGPROG =
endif
bin_PROGRAMS = lalapps_stringfind $(RINGPROG)
lalapps_stringfind_SOURCES = stringfind.c
if FRAME
if LALFRAME
if LALMETAIO
lalapps_ring_SOURCES = ring.c
TESTS = ring_test.sh
endif
endif
endif
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
MOSTLYCLEANFILES = ring-*.dat *RING*.xml *RING*.asc *RING*.gwf
EXTRA_DIST = ring.c stringfind.c \
             F-600000000-6-10.F F-600000060-6-10.F F-600000120-6-10.F \
             H-CAL_FAC-600000000-6000.gwf H-CAL_REF-600000000-64.gwf \
             H-600000000-600000180.cache H1-CAL-600000000-600006000.cache \
             response.asc filterpar.in ring_test.sh string_test.sh
