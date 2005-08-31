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
bin_SCRIPTS = lalapps_ring_pipe 
bin_PROGRAMS = $(RINGPROG)
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
endif
endif
endif
pythonlibdir = ${prefix}/lib/python
pythonlib_SCRIPTS = ring.py
pkgdata_DATA = ring_test.ini
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
CLEANFILES = python_config.sed lalapps_ring_pipe
EXTRA_DIST = ring.c $(SCRIPTS) $(DATA)

python_config.sed:
	-rm -f python_config.sed
	echo "s+@PYTHONLIBDIR@+${pythonlibdir}+g" > python_config.sed
	echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed
	echo "s+@RINGEXECPATH@+${bindir}+g" >> python_config.sed

lalapps_ring_pipe: ring_pipe.in python_config.sed
	sed -f python_config.sed ring_pipe.in > lalapps_ring_pipe 
	chmod +x lalapps_ring_pipe
