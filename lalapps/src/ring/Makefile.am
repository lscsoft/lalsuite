## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
if LALFRAME
if LALMETAIO
RINGPROG = lalapps_ring
STRINGPROG = 
else
RINGPROG =
STRINGPROG =
endif
else
RINGPROG =
STRINGPROG =
endif
else
RINGPROG =
STRINGPROG =
endif
bin_SCRIPTS = lalapps_ring_pipe
bin_PROGRAMS = $(RINGPROG) $(STRINGPROG)
if FRAME
if LALFRAME
if LALMETAIO
lalapps_ring_SOURCES = ring.c ring.h ring_option.c ring_filter.c ring_output.c \
         errutil.c errutil.h \
         gpstime.c gpstime.h \
	 getdata.c getdata.h \
	 injsgnl.c injsgnl.h \
	 getresp.c getresp.h \
	 spectrm.c spectrm.h \
	 segment.c segment.h
TESTS = ring_test.sh
lalapps_StringSearch_SOURCES = StringSearch.c
endif
endif
endif
pythonlibdir = ${prefix}/lib/python
pythonlib_SCRIPTS = ring.py
pkgdata_DATA = S2H1v03_test.txt S2_H1_dag.ini
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
MOSTLYCLEANFILES = ring-*.dat *RING*.xml *RING*.asc *RING*.gwf \
		   S2_H1_dag.*S2_RING_001.*
CLEANFILES = python_config.sed lalapps_ring_pipe S2_H1_dag.ini
EXTRA_DIST = ring.c StringSearch.c\
             F-600000000-6-10.F F-600000060-6-10.F F-600000120-6-10.F \
             H-CAL_FAC-600000000-6000.gwf H-CAL_REF-600000000-64.gwf \
             H-600000000-600000180.cache H1-CAL-600000000-600006000.cache \
             response.asc filterpar.in ring_test.sh \
	     ring_test_inj.sh ring_test_sim.sh ring_test_pipe.sh \
             $(SCRIPTS) $(DATA) \
	     ring_pipe.in S2_H1_dag.in S2H1v03_test.txt

EXTRA_PROGRAMS = lalapps_StringSearch

python_config.sed:
	-rm -f python_config.sed
	echo "s+@PYTHONLIBDIR@+${pythonlibdir}+g" > python_config.sed
	echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed
	echo "s+@RINGEXECPATH@+${bindir}+g" >> python_config.sed

lalapps_ring_pipe: ring_pipe.in python_config.sed
	sed -f python_config.sed ring_pipe.in > lalapps_ring_pipe 
	chmod +x lalapps_ring_pipe

S2_H1_dag.ini: S2_H1_dag.in python_config.sed
	sed "s+@RINGEXECPATH@+${bindir}+" S2_H1_dag.in > S2_H1_dag.ini 
