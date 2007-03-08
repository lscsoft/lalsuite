## Process this file with automake to produce Makefile.in
# $Id$
#
# Cristina
#
pythonlibdir = ${prefix}/lib/python
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps

TSUTILSRCS = tfdat2pgm.c tracksearchToolbox.c
TRACKSEARCHSRCS = tracksearch.c tracksearchToolbox.c
TSDATAGENSRCS = TSDatagen.c 
TRACKSEARCHAVGSRC = tracksearchAverager.c tracksearchToolbox.c

bin_PROGRAMS = lalapps_tracksearch lalapps_TSDatagen lalapps_TSaverager lalapps_tfdat2pgm
lalapps_tracksearch_SOURCES = $(TRACKSEARCHSRCS)
lalapps_TSDatagen_SOURCES = $(TSDATAGENSRCS)
lalapps_TSaverager_SOURCES = $(TRACKSEARCHAVGSRC)
lalapps_tfdat2pgm_SOURCES=$(TSUTILSRCS)

EXTRA_DIST = $(TRACKSEARCHSRCS)
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if HAVE_PYTHON
bin_SCRIPTS = candidateHandler.py tracksearch_pipe.py tuneTSpipe.py 
pkgpython_DATA = candidateUtils.py tracksearch.py
endif

pkgdata_DATA = tracksearchPipeline.ini tracksearchOpts.tun

