AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if LALFRAME
if LALBURST
bin_PROGRAMS = \
	lalapps_tracksearch \
	lalapps_tracksearch_datagen \
	lalapps_tracksearch_averager \
	lalapps_tracksearch_dat2pgm
endif
endif

noinst_HEADERS = \
	FrameDataConvert_source.c

lalapps_tracksearch_SOURCES = \
	tracksearch.c \
	tracksearchToolbox.c \
	tracksearch.h \
	tracksearchToolbox.h \
	tracksearchAverager.h \
	FrameDataConvert.c \
	FrameDataConvert.h

lalapps_tracksearch_datagen_SOURCES = \
	TSDatagen.c \
	TSDatagen.h \
	tracksearch.h \
	tracksearchToolbox.h \
	tracksearchAverager.h

lalapps_tracksearch_averager_SOURCES = \
	tracksearchAverager.c \
	tracksearchToolbox.c \
	tracksearch.h \
	tracksearchAverager.h \
	tracksearchToolbox.h

lalapps_tracksearch_dat2pgm_SOURCES = \
	tfdat2pgm.c \
	tracksearchToolbox.c \
	tracksearch.h \
	tracksearchAverager.h \
	tracksearchToolbox.h

if HAVE_PYTHON
bin_SCRIPTS = \
	lalapps_tracksearch_handler \
	lalapps_tracksearch_pipe \
	lalapps_tracksearch_housekeeper \
	lalapps_tracksearch_tune_pipe \
	lalapps_tracksearch_cutpipe \
	lalapps_tracksearch_mapcachebuilder
pkgpython_PYTHON = \
	tracksearchutils.py \
	tracksearch.py
endif

python_config.sed:
	@-rm -f python_config.sed
	@echo "s+@PYTHONLIBDIR@+@pkgpythondir@+g" > python_config.sed
	@echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed

lalapps_tracksearch_handler: $(srcdir)/tsHandler.pyin python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/tsHandler.pyin > lalapps_tracksearch_handler
	@chmod +x lalapps_tracksearch_handler

lalapps_tracksearch_pipe: $(srcdir)/tsPipe.pyin python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/tsPipe.pyin > lalapps_tracksearch_pipe
	@chmod +x lalapps_tracksearch_pipe

lalapps_tracksearch_cutpipe: $(srcdir)/tsCutPipe.pyin python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/tsCutPipe.pyin > lalapps_tracksearch_cutpipe
	@chmod +x lalapps_tracksearch_cutpipe

lalapps_tracksearch_housekeeper: $(srcdir)/tsHousekeeper.pyin python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/tsHousekeeper.pyin > lalapps_tracksearch_housekeeper
	@chmod +x lalapps_tracksearch_housekeeper

lalapps_tracksearch_tune_pipe: $(srcdir)/tsTunePipe.pyin python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/tsTunePipe.pyin > lalapps_tracksearch_tune_pipe
	@chmod +x lalapps_tracksearch_tune_pipe

lalapps_tracksearch_mapcachebuilder: $(srcdir)/tsMapCacheBuilder.pyin python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/tsMapCacheBuilder.pyin > lalapps_tracksearch_mapcachebuilder
	@chmod +x lalapps_tracksearch_mapcachebuilder

pkgdata_DATA = \
	tracksearch_pipeline_configuration.ini \
	tracksearch_tuning_configuration.tun \
	tracksearch_cutpipe_configuration.cut

EXTRA_DIST = \
	tsHandler.pyin \
	tsPipe.pyin \
	tsCutPipe.pyin \
	tsHousekeeper.pyin \
	tsTunePipe.pyin \
	tsMapCacheBuilder.pyin \
	$(DATA)

CLEANFILES = \
	lalapps_tracksearch_handler \
	lalapps_tracksearch_pipe \
	lalapps_tracksearch_housekeeper \
	lalapps_tracksearch_tune_pipe \
	lalapps_tracksearch_cutpipe \
	lalapps_tracksearch_mapcachebuilder \
	python_config.sed \
	lalapps_tracksearch \
	lalapps_tracksearch_datagen \
	lalapps_tracksearch_averager \
	lalapps_tracksearch_dat2pgmy

MOSTLYCLEANFILES = \
	*.candidates \
	*.dat \
	*.pgm \
	*.aux \
	*.raw \
	*.diag \
	*.pyc \
	tsaMap* \
	Candidate* \
	Pre-*

BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
