BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalapps.am

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
pybin_SCRIPTS = \
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

pkgdata_DATA = \
	tracksearch_pipeline_configuration.ini \
	tracksearch_tuning_configuration.tun \
	tracksearch_cutpipe_configuration.cut

EXTRA_DIST += \
	$(DATA)

CLEANFILES = \
	lalapps_tracksearch \
	lalapps_tracksearch_datagen \
	lalapps_tracksearch_averager \
	lalapps_tracksearch_dat2pgmy

MOSTLYCLEANFILES += \
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
