INCLUDES = -I$(top_srcdir)/src/lalapps -I$(top_srcdir)/misc

if LALFRAME
if LALMETAIO
bin_PROGRAMS = lalapps_power lalapps_binj
lalapps_power_SOURCES = power.c
lalapps_binj_SOURCES = binj.c
endif
endif

noinst_PROGRAMS = sense_compare
sense_compare_SOURCES = sense_compare.c

## make sure LALAPPSgitID.h is re-generated at every 'make'
.PHONY: gitID
gitID:
	cd $(top_builddir)/misc && $(MAKE) gitID
BUILT_SOURCES = gitID

if HAVE_PYTHON
bin_SCRIPTS = lalapps_power_pipe lalapps_power_online_pipe lalapps_power_likelihood_pipe
pkgpython_DATA = power.py
endif

pkgdata_DATA = power_pipe.ini power_online_pipe.ini

EXTRA_DIST = lalapps_power_pipe lalapps_power_online_pipe \
             lalapps_power_likelihood_pipe power.py \
             power_pipe.ini power_online_pipe.ini

LDADD = $(top_builddir)/src/lalapps/liblalapps.la
