## Process this file with automake to produce Makefile.in
# $Id$
POWERSRCS = power.c power.h 
BURCASRC = burca.c
BREADSRC = bread.c
BINJSRC = binj.c
BINJFINDSRC = binj_find.c

if LALFRAME
if LALMETAIO
bin_PROGRAMS = lalapps_power lalapps_burca \
               lalapps_bread lalapps_binj \
               lalapps_binj_find
lalapps_power_SOURCES = $(POWERSRCS)
lalapps_burca_SOURCES = $(BURCASRC)
lalapps_bread_SOURCES = $(BREADSRC)
lalapps_binj_SOURCES = $(BINJSRC)
lalapps_binj_find_SOURCES = $(BINJFINDSRC)
else
bin_PROGRAMS = 
endif
else
bin_PROGRAMS = 
endif
bin_SCRIPTS = lalapps_power_pipe 

pythonlibdir = ${prefix}/lib/python
pythonlib_SCRIPTS = power.py

lalapps_power_pipe: power_pipe.in
	sed s+@PYTHONLIBDIR@+${pythonlibdir}+g power_pipe.in > lalapps_power_pipe 

pkgdata_DATA = power_pipe.ini

EXTRA_DIST = $(POWERSRCS) power_pipe.in power_pipe.ini

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

CLEANFILES = lalapps_power_pipe 

MOSTLYCLEANFILES = *POWER*.xml *POWER*.gwf 
