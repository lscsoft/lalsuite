INCLUDES = -I$(top_srcdir)/src/lalapps

if LALFRAME
if LALMETAIO
bin_PROGRAMS = lalapps_power lalapps_bread lalapps_binj
lalapps_power_SOURCES = power.c
lalapps_bread_SOURCES = bread.c
lalapps_binj_SOURCES = binj.c
else
bin_PROGRAMS = 
endif
else
bin_PROGRAMS = 
endif

if HAVE_PYTHON
bin_SCRIPTS = lalapps_power_coinc_pipe lalapps_power_pipe lalapps_power_online_pipe
pkgpython_DATA = power.py
endif
EXTRA_DIST = power.c bread.c binj.c \
	     lalapps_power_coinc_pipe lalapps_power_pipe lalapps_power_online_pipe \
	     power.py power_pipe.ini power_online_pipe.ini

pkgdata_DATA = power_pipe.ini power_online_pipe.ini

LDADD = $(top_builddir)/src/lalapps/liblalapps.la
