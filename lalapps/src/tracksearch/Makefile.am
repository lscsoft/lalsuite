## Process this file with automake to produce Makefile.in
# $Id$
#
if LALFRAME
if LALMETAIO
bin_PROGRAMS = lalapps_tracksearch \
	lalapps_tracksearch_datagen \
	lalapps_tracksearch_averager \
	lalapps_tracksearch_dat2pgm
TSUTILSRCS = tfdat2pgm.c tracksearchToolbox.c tracksearch.h \
	tracksearchAverager.h tracksearchToolbox.h
TRACKSEARCHSRCS = tracksearch.c	tracksearchToolbox.c \
	tracksearch.h tracksearchToolbox.h tracksearchAverager.h
TSDATAGENSRCS = TSDatagen.c TSDatagen.h tracksearch.h \
	tracksearchToolbox.h tracksearchAverager.h
TRACKSEARCHAVGSRC = tracksearchAverager.c tracksearchToolbox.c \
	tracksearch.h tracksearchAverager.h tracksearchToolbox.h
lalapps_tracksearch_SOURCES = $(TRACKSEARCHSRCS)
lalapps_tracksearch_datagen_SOURCES = $(TSDATAGENSRCS)
lalapps_tracksearch_averager_SOURCES = $(TRACKSEARCHAVGSRC)
lalapps_tracksearch_dat2pgm_SOURCES = $(TSUTILSRCS)
else
bin_PROGRAMS =
endif
else
bin_PROGRAMS =
endif

if HAVE_PYTHON
bin_SCRIPTS = lalapps_tracksearch_handler \
	lalapps_tracksearch_pipe lalapps_tracksearch_housekeeper \
	lalapps_tracksearch_tune_pipe lalapps_tracksearch_cutpipe \
	lalapps_tracksearch_mapcachebuilder
pkgpython_DATA = tracksearchutils.py tracksearch.py
endif

python_config.sed:
	-rm -f python_config.sed
	echo "s+@PYTHONLIBDIR@+@pkgpythondir@+g" > python_config.sed
	echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed

lalapps_tracksearch_handler: $(srcdir)/tsHandler.pyin python_config.sed
	sed -f python_config.sed $(srcdir)/tsHandler.pyin > lalapps_tracksearch_handler
	chmod +x lalapps_tracksearch_handler

lalapps_tracksearch_pipe: $(srcdir)/tsPipe.pyin python_config.sed
	sed -f python_config.sed $(srcdir)/tsPipe.pyin > lalapps_tracksearch_pipe
	chmod +x lalapps_tracksearch_pipe

lalapps_tracksearch_cutpipe: $(srcdir)/tsCutPipe.pyin python_config.sed
	sed -f python_config.sed $(srcdir)/tsCutPipe.pyin > lalapps_tracksearch_cutpipe
	chmod +x lalapps_tracksearch_cutpipe

lalapps_tracksearch_housekeeper: $(srcdir)/tsHousekeeper.pyin python_config.sed
	sed -f python_config.sed $(srcdir)/tsHousekeeper.pyin > lalapps_tracksearch_housekeeper
	chmod +x lalapps_tracksearch_housekeeper

lalapps_tracksearch_tune_pipe: $(srcdir)/tsTunePipe.pyin python_config.sed
	sed -f python_config.sed $(srcdir)/tsTunePipe.pyin > lalapps_tracksearch_tune_pipe
	chmod +x lalapps_tracksearch_tune_pipe

lalapps_tracksearch_mapcachebuilder: $(srcdir)/tsMapCacheBuilder.pyin python_config.sed
	sed -f python_config.sed $(srcdir)/tsMapCacheBuilder.pyin > lalapps_tracksearch_mapcachebuilder
	chmod +x lalapps_tracksearch_mapcachebuilder

pkgdata_DATA = tracksearch_pipeline_configuration.ini \
	tracksearch_tuning_configuration.tun \
	tracksearch_cutpipe_configuration.cut
EXTRA_DIST = tsHandler.pyin tsPipe.pyin tsCutPipe.pyin \
	tsHousekeeper.pyin tsTunePipe.pyin tsMapCacheBuilder.pyin \
	$(DATA)

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

CLEANFILES = lalapps_tracksearch_handler \
	lalapps_tracksearch_pipe lalapps_tracksearch_housekeeper \
	lalapps_tracksearch_tune_pipe lalapps_tracksearch_cutpipe \
	lalapps_tracksearch_mapcachebuilder python_config.sed

MOSTLYCLEANFILES = *.candidates *.dat *.pgm *.aux *.raw *.diag \
	*.pyc tsaMap* Candidate* Pre-*
