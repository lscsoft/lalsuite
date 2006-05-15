## Process this file with automake to produce Makefile.in
# $Id$
POWERSRCS = power.c
BREADSRC = bread.c
BINJSRC = binj.c
COALESCENCEINFOSRC = coalescenceinfo.c

if LALFRAME
if LALMETAIO
bin_PROGRAMS = lalapps_power \
               lalapps_bread lalapps_binj \
               lalapps_coalescenceinfo
lalapps_power_SOURCES = $(POWERSRCS)
lalapps_bread_SOURCES = $(BREADSRC)
lalapps_binj_SOURCES = $(BINJSRC)
lalapps_coalescenceinfo_SOURCES = $(COALESCENCEINFOSRC)
else
bin_PROGRAMS = 
endif
else
bin_PROGRAMS = 
endif
bin_SCRIPTS = lalapps_power_coinc_pipe lalapps_power_pipe lalapps_power_online_pipe

pythonlibdir = ${prefix}/lib/python
pythonlib_SCRIPTS = power.py

pkgdata_DATA = power_pipe.ini

EXTRA_DIST = $(POWERSRCS) power_pipe.ini power.py

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

MOSTLYCLEANFILES = *POWER*.xml *POWER*.gwf 
