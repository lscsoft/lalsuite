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
pythonlibdir = ${prefix}/lib/python
pythonlib_SCRIPTS = ring.py
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
MOSTLYCLEANFILES = ring-*.dat *RING*.xml *RING*.asc *RING*.gwf
CLEANFILES = lalapps_ring_pipe
EXTRA_DIST = ring.c stringfind.c \
             F-600000000-6-10.F F-600000060-6-10.F F-600000120-6-10.F \
             H-CAL_FAC-600000000-6000.gwf H-CAL_REF-600000000-64.gwf \
             H-600000000-600000180.cache H1-CAL-600000000-600006000.cache \
             response.asc filterpar.in ring_test.sh string_test.sh \
             $(SCRIPTS) ring_pipe.in S2_H1_dag.ini

lalapps_ring_pipe: ring_pipe.in
	sed s+@PYTHONLIBDIR@+${pythonlibdir}+g ring_pipe.in > lalapps_ring_pipe 

