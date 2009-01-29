## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
if LALFRAME
if LALMETAIO
RINGPROG = lalapps_ring lalapps_rinj lalapps_ringread lalapps_rinca \
           lalapps_coincringread lalapps_collatesim
else
RINGPROG =
endif
else
RINGPROG =
endif
else
RINGPROG =
endif
bin_PROGRAMS = $(RINGPROG) \
	       lalapps_qnmode lalapps_qnmode_table lalapps_spheroidal_harmonic \
	       lalapps_ringdown_amplitude lalapps_ringdown_waveform
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
lalapps_rinj_SOURCES = rinj.c
lalapps_ringread_SOURCES = ringread.c
lalapps_rinca_SOURCES = rinca.c
lalapps_coincringread_SOURCES = coincringread.c
lalapps_collatesim_SOURCES = collatesim.c
endif
endif
endif
lalapps_qnmode_SOURCES = qnmode.c
lalapps_qnmode_table_SOURCES = qnmode_table.c
lalapps_spheroidal_harmonic_SOURCES = spheroidal_harmonic.c
lalapps_ringdown_amplitude_SOURCES = ringdown_amplitude.c
lalapps_ringdown_waveform_SOURCES = ringdown_waveform.c

## make sure LALAPPSgitID.h is re-generated at every 'make'
.PHONY: gitID
gitID:
	cd $(top_builddir)/misc && $(MAKE) gitID
BUILT_SOURCES = gitID

if HAVE_PYTHON
pkgpython_DATA = ring.py
bin_SCRIPTS = lalapps_ring_pipe 
endif

pkgdata_DATA = ring_test.ini
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_srcdir)/misc
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
CLEANFILES = python_config.sed lalapps_ring_pipe
EXTRA_DIST = ring.c ring_pipe.in ring.py $(DATA)

python_config.sed:
	-rm -f python_config.sed
	echo "s+@PYTHONLIBDIR@+@pkgpythondir@+g" > python_config.sed
	echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed
	echo "s+@RINGEXECPATH@+${bindir}+g" >> python_config.sed

lalapps_ring_pipe: $(srcdir)/ring_pipe.in python_config.sed
	sed -f python_config.sed $(srcdir)/ring_pipe.in > lalapps_ring_pipe 
	chmod +x lalapps_ring_pipe
