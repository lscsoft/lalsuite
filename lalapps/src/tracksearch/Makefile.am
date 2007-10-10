## Process this file with automake to produce Makefile.in
# $Id$
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

EXTRA_DIST = $(TRACKSEARCHSRCS)
pkgdata_DATA = tracksearchPipeline.ini tracksearchOpts.tun
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if HAVE_PYTHON
bin_SCRIPTS = lalapps_tracksearch_handler \
	      lalapps_tracksearch_pipe \
              lalapps_tracksearch_housekeeper \
              lalapps_tracksearch_tune_pipe
pkgpython_DATA = tracksearchutils.py tracksearch.py
endif

python_config.sed:
	-rm -f python_config.sed
	echo "s+@PYTHONLIBDIR@+@pkgpythondir@+g" > python_config.sed
	echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed

lalapps_tracksearch_handler: $(srcdir)/tshandler python_config.sed
	sed -f python_config.sed $(srcdir)/tshandler > lalapps_tracksearch_handler
	chmod +x lalapps_tracksearch_handler

lalapps_tracksearch_pipe: $(srcdir)/tspipe python_config.sed
	sed -f python_config.sed $(srcdir)/tspipe > lalapps_tracksearch_pipe
	chmod +x lalapps_tracksearch_pipe

lalapps_tracksearch_housekeeper: $(srcdir)/tspipe python_config.sed
	sed -f python_config.sed $(srcdir)/tspipe > lalapps_tracksearch_housekeeper
	chmod +x lalapps_tracksearch_housekeeper

lalapps_tracksearch_tune_pipe: $(srcdir)/tstunepipe python_config.sed
	sed -f python_config.sed $(srcdir)/tstunepipe > lalapps_tracksearch_tune_pipe
	chmod +x lalapps_tracksearch_tune_pipe

CLEANFILES =  lalapps_tracksearch_handler \
	      lalapps_tracksearch_pipe \
              lalapps_tracksearch_housekeeper \
              lalapps_tracksearch_tune_pipe

MOSTLYCLEANFILES = *.candidates \
                   *.dat \
                   *.pgm \
                   *.raw \
                   *.diag \
                   tsaMap* \
                   Candidate* \
                   Pre-*


