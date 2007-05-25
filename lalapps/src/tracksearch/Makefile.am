## Process this file with automake to produce Makefile.in
# $Id$
#
# Cristina
#
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
if LALFRAME
if LALMETAIO
bin_PROGRAMS = lalapps_tracksearch lalapps_TSDatagen lalapps_TSaverager lalapps_tfdat2pgm
TSUTILSRCS = tfdat2pgm.c tracksearchToolbox.c
TRACKSEARCHSRCS = tracksearch.c tracksearchToolbox.c
TSDATAGENSRCS = TSDatagen.c 
TRACKSEARCHAVGSRC = tracksearchAverager.c tracksearchToolbox.c
lalapps_tracksearch_SOURCES = $(TRACKSEARCHSRCS)
lalapps_TSDatagen_SOURCES = $(TSDATAGENSRCS)
lalapps_TSaverager_SOURCES = $(TRACKSEARCHAVGSRC)
lalapps_tfdat2pgm_SOURCES=$(TSUTILSRCS)
else
bin_PROGRAMS =
endif
else
bin_PROGRAMS =
endif

if HAVE_PYTHON
bin_SCRIPTS = candidateHandler.py tracksearch_pipe.py tuneTSpipe.py 
pkgpython_DATA = candidateUtils.py tracksearch.py
endif

EXTRA_DIST = $(TRACKSEARCHSRCS)
pkgdata_DATA = tracksearchPipeline.ini tracksearchOpts.tun
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
